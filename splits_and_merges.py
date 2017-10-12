import pygplates
import numpy as np


def get_topologies_and_plate_id_list(topology_features,rotation_model,time):
    resolved_topologies = []
    pygplates.resolve_topologies(topology_features, rotation_model, resolved_topologies, time)

    plate_id_list = []
    for topology in resolved_topologies:
        plate_id_list.append(topology.get_feature().get_reconstruction_plate_id())
        
    plate_area_list = []
    for topology in resolved_topologies:
        plate_area_list.append(topology.get_resolved_geometry().get_area())
        
    return resolved_topologies,plate_id_list,plate_area_list


def find_plate_change_for_time_step(topology_features,rotation_model,time_from,time_to,verbose=True):

    (resolved_topologies_t1,
     plate_id_list_t1,
     plate_area_list_t1) = get_topologies_and_plate_id_list(topology_features,rotation_model,time_to)
    (resolved_topologies_t2,
     plate_id_list_t2,
     plate_area_list_t2) = get_topologies_and_plate_id_list(topology_features,rotation_model,time_from)

    # compare lists from two time snapshots
    common_plates = set(plate_id_list_t1).intersection(plate_id_list_t2)
    appearing_plates = set(plate_id_list_t1).difference(plate_id_list_t2)
    disappearing_plates = set(plate_id_list_t2).difference(plate_id_list_t1)

    if verbose:
        # plate_ids in both lists
        print 'plates that persist between %0.2f Ma and %0.2f Ma are: \n %s\n' % \
            (time_from,time_to,common_plates)

        # plates ids that are not in the t2 list
        print 'plates that appeared between %0.2f Ma and %0.2f Ma are: \n %s\n' % \
            (time_from,time_to,appearing_plates)

        # plate ids that are not in the t1 list
        print 'plates that disappeared between %0.2f Ma and %0.2f Ma are: \n %s' % \
            (time_from,time_to,disappearing_plates)
        
        print '\nWorking on appearing plates'
    appearing_plate_map = plate_change_mapping(appearing_plates,
                                               resolved_topologies_t1,
                                               resolved_topologies_t2,
                                               verbose=verbose)
    
    if verbose:
        print '\nWorking on disappearing plates'
    disappearing_plate_map = plate_change_mapping(disappearing_plates,
                                                  resolved_topologies_t2,
                                                  resolved_topologies_t1,
                                                  verbose=verbose)
    #print '\n'
    
    return appearing_plate_map,disappearing_plate_map


def plate_change_mapping(plate_list_at_time,topologies_at_time,topologies_at_delta_time,verbose=True):
    
    plate_map = []
    
    for plate in plate_list_at_time:
        for topology in topologies_at_time:
            if topology.get_feature().get_reconstruction_plate_id()==plate:
                centroid = topology.get_resolved_geometry().get_interior_centroid()

        for topology in topologies_at_delta_time:
            if topology.get_resolved_geometry().is_point_in_polygon(centroid):
                if verbose:
                    print 'Centroid for plate %d mapped to plate %d at delta time' % \
                        (plate,topology.get_feature().get_reconstruction_plate_id())    
                plate_map.append((plate,topology.get_feature().get_reconstruction_plate_id()))
                    
    return plate_map



#def determine(moving_plate_id,fixed_plate_id,time_from,time_to,rotation_model):
#    angle = rotation_model.get_rotation(time_to,moving_plate_id,time_from,fixed_plate_id)
#    return angle.represents_identity_rotation()
#plate_id_list_copy_True[:] = [tup for tup in plate_id_list_copy_True if determine(tup,fixed_plate_id,time_from,time_to,rotation_model)]

def match_splits_to_origin(appearing_plate_map):
    # this function only cares about appearing plates, and aims
    # to find pairs of plates that split from the same common plate

    origin_plates_for_new_plates = []
    for plates in appearing_plate_map:
        origin_plates_for_new_plates.append(plates[1])


    new_plates_grouped_by_origin_plate = []
    for origin_plate in set(origin_plates_for_new_plates):
        successor_plates = []
        for plates in appearing_plate_map:
            if plates[1]==origin_plate:
                successor_plates.append(plates[0])
        new_plates_grouped_by_origin_plate.append((origin_plate,successor_plates))

    return new_plates_grouped_by_origin_plate
  
    

# This assumes that the plates always split into two (and not three, four.....)
def get_great_circle_along_plate_split(plate_split,
                                       shared_boundary_sections):
    
    great_circle_arc = None

    print plate_split
    geometries_along_new_split = []

    for shared_boundary_section in shared_boundary_sections:
        for shared_sub_segment in shared_boundary_section.get_shared_sub_segments():
            sharing_topologies = shared_sub_segment.get_sharing_resolved_topologies()

            if len(sharing_topologies)!=2:
                #print 'bad topological segment'
                continue
            else:    
                sub_segment_plate_pair = []
                for topology in sharing_topologies:
                    sub_segment_plate_pair.append(topology.get_feature().get_reconstruction_plate_id())
                #print sub_segment_plate_pair
                if (plate_split[1][0] in sub_segment_plate_pair and plate_split[1][1] in sub_segment_plate_pair):
                    geometries_along_new_split.append(shared_sub_segment.get_geometry())

    if len(geometries_along_new_split)!=1:
        print 'Only works for new splits with one geometry so far'
    else:
        #print geometries_along_new_split
        great_circle_arc = pygplates.GreatCircleArc(
            geometries_along_new_split[0].get_points()[0],
            geometries_along_new_split[0].get_points()[-1])  
                
    return great_circle_arc


def get_plate_disappearance_time_lut(topology_features,
                                     rotation_model,
                                     time_list,
                                     verbose=False):

    plate_disappearance_time_lut = []

    for time in time_list:

        time_from = time+1
        time_to = time

        (appearing_plate_map,
         disappearing_plate_map) = find_plate_change_for_time_step(topology_features,rotation_model,time_from,time_to,verbose)

        if len(zip(*disappearing_plate_map))>0:
            print list(zip(*disappearing_plate_map)[0])
            for plate_id in list(zip(*disappearing_plate_map)[0]):
                plate_id
                plate_disappearance_time_lut.append((plate_id,time))

    return plate_disappearance_time_lut

