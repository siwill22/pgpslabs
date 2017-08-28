from __future__ import print_function
import pygplates
import numpy as np
import sys
import xarray as xr
import inpaint
import scipy.interpolate as spi


bug_fix_arc_dir = -1.0 if pygplates.Version.get_imported_version() < pygplates.Version(14) else 1.0


# function that takes a netcdf grid, fills dummy values, then creates 
# an interpolator object that can be evaluated later at specified points
def make_age_interpolator(grdfile,interp='Spherical'):
 
    ds_disk = xr.open_dataset(grdfile)
    
    data_array = ds_disk['z']

    coord_keys = data_array.coords.keys()
    gridX = data_array.coords[coord_keys[0]].data
    gridY = data_array.coords[coord_keys[1]].data
    gridZ = data_array.data
        
    gridZ_filled = inpaint.fill_ndimage(gridZ)
    
    # spherical interpolation
    # Note the not-ideal method for avoiding issues with points at the edges of the grid
    if interp is 'Spherical':
        lut = spi.RectSphereBivariateSpline(np.radians(gridY[1:-1]+90.),
                                            np.radians(gridX[1:-1]+180.),
                                            gridZ_filled[1:-1,1:-1])

    # flat earth interpolation
    elif interp is 'FlatEarth':
        lut=spi.RectBivariateSpline(gridX,gridY,gridZ_filled.T)

    return lut


def getSubductionBoundarySections(topology_features,rotation_model,time):
# given files to make topological polygons, returns the features of type 'MidOceanRidge'
# and get the first and last point from each one, along with the plate pairs
    
    subduction_boundary_sections = []
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time, shared_boundary_sections)
                    
    for shared_boundary_section in shared_boundary_sections:
                
        if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('SubductionZone'):
                 
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                
                subduction_boundary_sections.append(shared_sub_segment)
                              
    return subduction_boundary_sections


# Determine the overriding and subducting plates of the subduction shared sub-segment.
def find_overriding_and_subducting_plates(subduction_shared_sub_segment, time):
    
    # Get the subduction polarity of the nearest subducting line.
    subduction_polarity = subduction_shared_sub_segment.get_feature().get_enumeration(pygplates.PropertyName.gpml_subduction_polarity)
    if (not subduction_polarity) or (subduction_polarity == 'Unknown'):
        print('Unable to find the overriding plate of the subducting shared sub-segment "{0}"'.format(
            subduction_shared_sub_segment.get_feature().get_name()), file=sys.stderr)
        print('    subduction zone feature is missing subduction polarity property or it is set to "Unknown".', file=sys.stderr)
        return

    # There should be two sharing topologies - one is the overriding plate and the other the subducting plate.
    sharing_resolved_topologies = subduction_shared_sub_segment.get_sharing_resolved_topologies()
    if len(sharing_resolved_topologies) != 2:
        print('Unable to find the overriding and subducting plates of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    there are not exactly 2 topologies sharing the sub-segment.', file=sys.stderr)
        return

    overriding_plate = None
    subducting_plate = None
    
    geometry_reversal_flags = subduction_shared_sub_segment.get_sharing_resolved_topology_geometry_reversal_flags()
    for index in range(2):

        sharing_resolved_topology = sharing_resolved_topologies[index]
        geometry_reversal_flag = geometry_reversal_flags[index]

        if sharing_resolved_topology.get_resolved_boundary().get_orientation() == pygplates.PolygonOnSphere.Orientation.clockwise:
            # The current topology sharing the subducting line has clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line not reversed in the topology.
            if ((subduction_polarity == 'Left' and geometry_reversal_flag) or
                (subduction_polarity == 'Right' and not geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
        else:
            # The current topology sharing the subducting line has counter-clockwise orientation (when viewed from above the Earth).
            # If the overriding plate is to the 'left' of the subducting line (when following its vertices in order) and
            # the subducting line is not reversed when contributing to the topology then that topology is the overriding plate.
            # A similar test applies to the 'right' but with the subducting line reversed in the topology.
            if ((subduction_polarity == 'Left' and not geometry_reversal_flag) or
                (subduction_polarity == 'Right' and geometry_reversal_flag)):
                overriding_plate = sharing_resolved_topology
            else:
                subducting_plate = sharing_resolved_topology
    
    if overriding_plate is None:
        print('Unable to find the overriding plate of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    both sharing topologies are on subducting side of subducting line.', file=sys.stderr)
        return
    
    if subducting_plate is None:
        print('Unable to find the subducting plate of the subducting shared sub-segment "{0}" at {1}Ma'.format(
            subduction_shared_sub_segment.get_feature().get_name(), time), file=sys.stderr)
        print('    both sharing topologies are on overriding side of subducting line.', file=sys.stderr)
        return
    
    return (overriding_plate, subducting_plate, subduction_polarity)


# function to warp polyline based on starting subduction segment location
def warp_subduction_segment(tessellated_line,
                            rotation_model,
                            subducting_plate_id,
                            subduction_polarity,
                            time,
                            end_time,
                            time_step,
                            dip_angle_radians,
                            use_small_circle_path=False):

    # We need to reverse the subducting_normal vector direction if overriding plate is to
    # the right of the subducting line since great circle arc normal is always to the left.
    if subduction_polarity == 'Left':
        subducting_normal_reversal = 1
    else:
        subducting_normal_reversal = -1

    # tesselate the line, and create an array with the same length as the tesselated points
    # with zero as the starting depth for each point at t=0
    points = [point for point in tessellated_line]
    point_depths = [0. for point in points]

    # Need at least two points for a polyline. Otherwise, return None for all results
    if len(points) < 2:
        points = None; point_depths = None; polyline = None
        return points, point_depths, polyline

    polyline = pygplates.PolylineOnSphere(points)

    warped_polylines = []

    # Add original unwarped polyline first.
    warped_polylines.append(polyline)

    #warped_end_time = time - warped_time_interval
    warped_end_time = end_time
    if warped_end_time < 0:
        warped_end_time = 0

    # iterate over each time in the range defined by the input parameters
    for warped_time in np.arange(time, warped_end_time-time_step,-time_step):

        # the stage rotation that describes the motion of the subducting plate,
        # with respect to the fixed plate for the rotation model
        stage_rotation = rotation_model.get_rotation(warped_time-time_step, subducting_plate_id, warped_time)

        if use_small_circle_path:
            stage_pole, stage_pole_angle_radians = stage_rotation.get_euler_pole_and_angle()

        # get velocity vectors at each point along polyline
        relative_velocity_vectors = pygplates.calculate_velocities(
                points,
                stage_rotation,
                time_step,
                pygplates.VelocityUnits.kms_per_my)

        # Get subducting normals for each segment of tessellated polyline.
        # Also add an imaginary normal prior to first and post last points
        # (makes its easier to later calculate average normal at tessellated points).
        # The number of normals will be one greater than the number of points.
        subducting_normals = []
        subducting_normals.append(None) # Imaginary segment prior to first point.
        for segment in polyline.get_segments():
            if segment.is_zero_length():
                subducting_normals.append(None)
            else:
                # The normal to the subduction zone in the direction of subduction (towards overriding plate).
                subducting_normals.append(
                    subducting_normal_reversal * segment.get_great_circle_normal())
        subducting_normals.append(None) # Imaginary segment after to last point.

        # get vectors of normals and parallels for each segment, use these 
        # to get a normal and parallel at each point location
        normals = []
        parallels = []            
        for point_index in range(len(points)):
            prev_normal = subducting_normals[point_index]
            next_normal = subducting_normals[point_index + 1]

            if prev_normal is None and next_normal is None:
                # Skip point altogether (both adjoining segments are zero length).
                continue

            if prev_normal is None:
                normal = next_normal
            elif next_normal is None:
                normal = prev_normal
            else:
                normal = (prev_normal + next_normal).to_normalised()

            parallel = pygplates.Vector3D.cross(point.to_xyz(), normal).to_normalised()

            normals.append(normal)
            parallels.append(parallel)

        # iterate over each point to determine the incremented position 
        # based on plate motion and subduction dip
        warped_points = []
        warped_point_depths = []
        for point_index, point in enumerate(points):
            normal = normals[point_index]
            parallel = parallels[point_index]

            velocity = relative_velocity_vectors[point_index]
            if velocity.is_zero_magnitude():
                # Point hasn't moved.
                warped_points.append(point)
                warped_point_depths.append(point_depths[point_index])
                continue

            # reconstruct the tracked point from position at current time to
            # position at the next time step
            normal_angle = pygplates.Vector3D.angle_between(velocity, normal)
            parallel_angle = pygplates.Vector3D.angle_between(velocity, parallel)

            # Trench parallel and normal components of velocity.
            velocity_normal = np.cos(normal_angle) * velocity.get_magnitude()
            velocity_parallel = np.cos(parallel_angle) * velocity.get_magnitude()

            normal_vector = normal.to_normalised() * velocity_normal
            parallel_vector = parallel.to_normalised() * velocity_parallel

            # Adjust velocity based on subduction vertical dip angle.
            velocity_dip = parallel_vector + np.cos(dip_angle_radians) * normal_vector

            #deltaZ is the amount that this point increases in depth within the time step
            deltaZ = np.sin(dip_angle_radians) * velocity.get_magnitude()

            # Should be 90 degrees always.
            #print np.degrees(np.arccos(pygplates.Vector3D.dot(normal_vector, parallel_vector)))

            if use_small_circle_path:
                # Rotate original stage pole by the same angle that effectively
                # rotates the velocity vector to the dip velocity vector.
                dip_stage_pole_rotate = pygplates.FiniteRotation(
                        point,
                        pygplates.Vector3D.angle_between(velocity_dip, velocity))
                dip_stage_pole = dip_stage_pole_rotate * stage_pole
            else:
                # Get the unnormalised vector perpendicular to both the point and velocity vector.
                dip_stage_pole_x, dip_stage_pole_y, dip_stage_pole_z = pygplates.Vector3D.cross(
                        point.to_xyz(), velocity_dip).to_xyz()

                # PointOnSphere requires a normalised (ie, unit length) vector (x, y, z).
                dip_stage_pole = pygplates.PointOnSphere(
                        dip_stage_pole_x, dip_stage_pole_y, dip_stage_pole_z, normalise=True)


            # Get angle that velocity will rotate seed point along great circle arc 
            # over 'time_step' My (if velocity in Kms / My).
            dip_stage_angle_radians = velocity_dip.get_magnitude() * (
                    time_step / pygplates.Earth.mean_radius_in_kms)

            if use_small_circle_path:
                # Increase rotation angle to adjust for fact that we're moving a
                # shorter distance with small circle (compared to great circle).
                dip_stage_angle_radians /= np.abs(np.sin(
                        pygplates.Vector3D.angle_between(
                                dip_stage_pole.to_xyz(), point.to_xyz())))
                # Use same sign as original stage rotation.
                if stage_pole_angle_radians < 0:
                    dip_stage_angle_radians = -dip_stage_angle_radians

            # get the stage rotation that describes the lateral motion of the 
            # point taking the dip into account
            dip_stage_rotation = pygplates.FiniteRotation(dip_stage_pole, dip_stage_angle_radians)

            # increment the point long,lat and depth
            warped_point = dip_stage_rotation * point
            warped_points.append(warped_point)
            warped_point_depths.append(point_depths[point_index] + deltaZ)

        # finished warping all points in polyline
        # --> increment the polyline for this time step
        warped_polyline = pygplates.PolylineOnSphere(warped_points)
        warped_polylines.append(warped_polyline)

        # For next warping iteration.
        points = warped_points
        polyline = warped_polyline
        point_depths = warped_point_depths
        
    return points, point_depths, polyline


def write_subducted_slabs_to_xyz(output_filename,output_data):
	
    with open(output_filename, 'w') as output_file:
        output_file.write('Long,Lat,Depth,AgeAtSubduction,TimeOfSubduction\n')
        for output_segment in output_data:
            for index in range(len(output_segment[2])):
                output_file.write('%0.6f,%0.6f,%0.6f,%0.2f,%0.2f\n' % (output_segment[1].to_lat_lon_array()[index,1],
                                                                       output_segment[1].to_lat_lon_array()[index,0],
                                                                       output_segment[2][index],
                                                                       output_segment[3][index],
                                                                       output_segment[0]))



#######################################################################
# Function from here downwards are largely deprecated
def getRidgeEndPoints(topology_features,rotation_model,time):
# given files to make topological polygons, returns the features of type 'MidOceanRidge'
# and get the first and last point from each one, along with the plate pairs
    
    MorEndPointArrays = []
    MorEndPointGeometries = []
    MorPlatePairs = []
    subduction_boundary_sections = []
    
    # Resolve our topological plate polygons (and deforming networks) to the current 'time'.
    # We generate both the resolved topology boundaries and the boundary sections between them.
    resolved_topologies = []
    shared_boundary_sections = []
    pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time, shared_boundary_sections)
                    
    for shared_boundary_section in shared_boundary_sections:
        if shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('MidOceanRidge'):
                 
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                
                if len(shared_sub_segment.get_sharing_resolved_topologies())==2:
                    plate_pair = [shared_sub_segment.get_sharing_resolved_topologies()[0].get_feature().get_reconstruction_plate_id(),
                                  shared_sub_segment.get_sharing_resolved_topologies()[1].get_feature().get_reconstruction_plate_id()]
                else:
                    plate_pair = np.array((-1,-1))
                #print 'skipping bad topological segment....'

                tmp = shared_sub_segment.get_geometry()
                
                MorEndPointArrays.append(tmp.to_lat_lon_array())
                MorEndPointGeometries.append(pygplates.PointOnSphere(tmp.get_points()[0]))
                MorEndPointGeometries.append(pygplates.PointOnSphere(tmp.get_points()[-1]))
                MorPlatePairs.append(plate_pair)
                MorPlatePairs.append(plate_pair)
                
        elif shared_boundary_section.get_feature().get_feature_type() == pygplates.FeatureType.create_gpml('SubductionZone'):
                 
            for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
                
                subduction_boundary_sections.append(shared_sub_segment)
                
                
    return MorEndPointArrays,MorEndPointGeometries,MorPlatePairs,subduction_boundary_sections


def getMor2szDistance2(MorEndPointGeometries,MorPlatePairs,subduction_boundary_sections):
# Get distance between end points of resolved mid-ocean ridge feature lines
# and subduction zones, using subduction segments from resolve_topologies

    Mor2szDistance = []
    sz_opid = []
    
    for MorPoint,PlatePair in zip(MorEndPointGeometries,MorPlatePairs):

        min_distance_to_all_features = np.radians(180)
        #nearest_sz_point = None
        opid = 0

        for subduction_boundary_section in subduction_boundary_sections:

            if MorPoint is not None:
                min_distance_to_feature = pygplates.GeometryOnSphere.distance(
                    MorPoint,
                    subduction_boundary_section.get_geometry(),
                    min_distance_to_all_features)

                # If the current geometry is nearer than all previous geometries then
                # its associated feature is the nearest feature so far.
                if min_distance_to_feature is not None:
                    min_distance_to_all_features = min_distance_to_feature
                    opid = subduction_boundary_section.get_feature().get_reconstruction_plate_id()

        Mor2szDistance.append(min_distance_to_all_features*pygplates.Earth.mean_radius_in_kms)
        sz_opid.append(opid)
        
    return Mor2szDistance,sz_opid

def track_point_to_present_day(seed_geometry,PlateID,rotation_model,start_time,end_time,time_step):
# Given a seed geometry at some time in the past, return the locations of this point
# at a series of times between that time and present-day

    point_longitude = []
    point_latitude = []
    
    for time in np.arange(start_time,end_time,-time_step):
        
        stage_rotation = rotation_model.get_rotation(time-time_step, PlateID, time)

        # use the stage rotation to reconstruct the tracked point from position at current time 
        # to position at the next time step
        incremented_geometry = stage_rotation * seed_geometry

        # replace the seed point geometry with the incremented geometry in preparation for next iteration
        seed_geometry = incremented_geometry

        point_longitude.append(seed_geometry.to_lat_lon_point().get_longitude())
        point_latitude.append(seed_geometry.to_lat_lon_point().get_latitude())
        
    return point_longitude,point_latitude


def rotate_point_to_present_day(seed_geometry,PlateID,rotation_model,start_time):
# Given a seed geometry at some time in the past, return the locations of this point
# at present-day (only)
   
    point_longitude = []
    point_latitude = []
    
    stage_rotation = rotation_model.get_rotation(0, PlateID, time, anchor_plate_id=1)

    # use the stage rotation to reconstruct the tracked point from position at current time 
    # to position at the next time step
    incremented_geometry = stage_rotation * seed_geometry

    # replace the seed point geometry with the incremented geometry in preparation for next iteration
    seed_geometry = incremented_geometry

    point_longitude.append(seed_geometry.to_lat_lon_point().get_longitude())
    point_latitude.append(seed_geometry.to_lat_lon_point().get_latitude())
        
    return point_longitude,point_latitude

def create_seed_point_feature(plat,plon,plate_id,conjugate_plate_id,time):
# Create a gpml point feature given some attributes 

    point = pygplates.PointOnSphere(plat,plon)
    point_feature = pygplates.Feature()
    point_feature.set_geometry(point)
    point_feature.set_valid_time(time,0.)
    point_feature.set_reconstruction_plate_id(plate_id)
    point_feature.set_conjugate_plate_id(conjugate_plate_id)
    point_feature.set_name('Slab Edge | plate %d | rel. plate %d' % (plate_id,conjugate_plate_id))
    
    return point_feature

