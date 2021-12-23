###########################################################
#
# Torino subtraction scheme -- Dipole Current class
#
###########################################################

import commons.general_current as general_current
import madgraph.integrator.mappings as mappings
import madgraph.core.subtraction as sub
import bidict
import torino_config as TRN_config
import factors_and_cuts as fc


import logging
logger1 = logging.getLogger('madgraph')


class DipoleCurrent(general_current.GeneralCurrent):
    """A new daughter class which makes it possible to
    specify different recoilers for different color dipoles
    """

    def get_mapping_info_for_all_steps(self, reduced_process, global_variables, *args, **opts):
        """ This function returns a list of tuples of the following format:
             ( reduced_kinematic_identifier, (
                                        mapping_info_for_step_1,
                                        mapping_info_for_step_2,
                                        ...
                                        )
            )
            By default this returns ( None, []) and the list of mapping information from the entries 'reduced_recoilers',
            'mapping_class' and 'additional_recoilers' of the mapping rules.
            Note that the leg numbers of these set of recoilers specified should be already mapped, i.e correspond
            to the leg numbers of the actual PS point and process supplied at the time this current is called.

            The mapping_info specified above are dictionaries of the following form:

                {
                    'mapping_class' : <class_that_implements_the_mapping>,
                    'mapping_singular_structure' : <singular_structure_specifying_the_input_to_the_mapping_and_recoilers>,
                    'mapping_momenta_dict: <bidictionary_indicating_the_momenta_label_routing_that_the_mapping_must_adopt>
                }

            For Catani-Seymour dipole subtraction for example, you would use a different mapping for each dipole, so one
            could return a list of entries like the following for the dipole with legs (1,2):
            (   (1,2),   [
                    'mapping_class': mapping.FinalLorentzOneMapping,
                    'mapping_singular_structure': sub.SingularStrcuture(
                        substructure=sub.CollStructure(legs=SubtractionLegSet(SubtractionLeg(number=1))),
                        legs=SubtractionLegSet(SubtractionLeg(number=2)
                    ),
                    'mapping_momenta_dict' : bidict({frozenset(1,2):3})
                ]
            )
        """

        mapping_info_list = []

        sector_legs = global_variables['sector_info'][0].leg_numbers

        overall_parents = global_variables['overall_parents']
        leg_numbers_map = global_variables['leg_numbers_map']

#        logger1.info('sector_legs : ' + str(sector_legs))
#        logger1.info('overall_parents : ' + str(overall_parents))
#        logger1.info('leg_numbers_map : ' + str(leg_numbers_map))

        #logger1.info('reduced_process : ' + str(global_variables))

        for i_step, mapping_information in enumerate(self.mapping_rules):

            # Now recursively apply leg numbers mappings
            mapping_singular_structure = mapping_information['singular_structure'].get_copy()
            self.map_leg_numbers_in_singular_structure(leg_numbers_map, mapping_singular_structure, overall_parents)
            singular_structure = mapping_singular_structure.substructures[0]
            singular_legs = singular_structure.legs

            for ia, lega in self.get_colored_legs(reduced_process,global_variables).items():
                for ib, legb in self.get_colored_legs(reduced_process,global_variables).items():
                    if ia == ib:
                        continue

                    # identify the leg that belongs to the sector
                    # and the one that does not
                    sec_leg = lega
                    oth_leg = legb
#                    if ia in sector_legs:
#                        sec_leg = lega
#                        oth_leg = legb
#                    elif ib in sector_legs:
#                        sec_leg = legb
#                        oth_leg = lega
#                    else:
#                        # in this case, the leg assignment is arbitrary.
#                        sec_leg = lega
#                        oth_leg = legb

#                    logger1.info('singular_structure name : ' + str(singular_structure.name()))
                    #logger1.info('ia : ' + str(lega))
                    #logger1.info('ib : ' + str(legb))


#                    if singular_structure.name() == 'C' and singular_structure.substructures[0].name() == 'S':
#                        logger1.info('Here we are: ' + str(singular_structure))

#                        mapping_structure = sub.SingularStructure(substructures=(singular_structure, ), \
#                                                        legs=sub.SubtractionLegSet((oth_leg,)))
#                        logger1.info('Legs : ' + str(sec_leg) + ',' + str(oth_leg)) 
#                        logger1.info('self.map_leg_number(leg_numbers_map, k, overall_parents) : ' + str(mapping_information['momenta_dict'].items())) 

#                        mapping_dict = bidict({self.map_leg_number(leg_numbers_map, k, overall_parents): frozenset([
#                            self.map_leg_number(leg_numbers_map, n, overall_parents) for n in v]) for k, v in
#                            mapping_information['momenta_dict'].items()})


#                    elif singular_structure.name() == 'C' and singular_structure.substructures[0].name() != 'S':
                    if singular_structure.name() == 'C' and ia == 1 and ib == 2:
                        # collinear structure: one among ia and ib is already
                        # in the singular_structure. Add the other
                        mapping_structure = sub.SingularStructure(substructures=(singular_structure, ), \
                                                        legs=sub.SubtractionLegSet((oth_leg,)))

                        # Build the momenta_dict by also substituting leg numbers
#                        mapping_dict = bidict.bidict({self.map_leg_number(leg_numbers_map, k, overall_parents): frozenset([
#                            self.map_leg_number(leg_numbers_map, 11, overall_parents), ia]) for k, v in
#                            mapping_information['momenta_dict'].items()})
                     

                        mapping_dict = bidict.bidict({self.map_leg_number(leg_numbers_map, k, overall_parents): frozenset([
                            self.map_leg_number(leg_numbers_map, n, overall_parents) for n in v]) for k, v in
                            mapping_information['momenta_dict'].items()})

                    elif singular_structure.name() == 'C' and ia == 2 and ib == 1:
                        # collinear structure: one among ia and ib is already
                        # in the singular_structure. Add the other
                        mapping_structure = sub.SingularStructure(substructures=(singular_structure, ), \
                                                        legs=sub.SubtractionLegSet((oth_leg,)))

                        # Build the momenta_dict by also substituting leg numbers
                        mapping_dict = bidict.bidict({self.map_leg_number(leg_numbers_map, k, overall_parents): frozenset([
                            self.map_leg_number(leg_numbers_map, n, overall_parents) for n in v]) for k, v in
                            mapping_information['momenta_dict'].items()})

                    elif singular_structure.name() == 'C':
                        #logger1.info('Extra collinear in dipole current')
                        # collinear structure: one among ia and ib is already
                        # in the singular_structure. Add the other
                        mapping_structure = sub.SingularStructure(substructures=(singular_structure, ), \
                                                        legs=sub.SubtractionLegSet((oth_leg,)))

                        # Build the momenta_dict by also substituting leg numbers
                        mapping_dict = bidict.bidict({self.map_leg_number(leg_numbers_map, k, overall_parents): frozenset([
                            self.map_leg_number(leg_numbers_map, n, overall_parents) for n in v]) for k, v in
                            mapping_information['momenta_dict'].items()})


                    elif singular_structure.name() == 'S':
                        # soft structure. Build first a collinear structure
                        # with sec_leg, then add oth_leg
                        mapping_structure = sub.SingularStructure(substructures=(
                                              sub.CollStructure(substructures=(singular_structure, ), \
                                                                legs=sub.SubtractionLegSet((sec_leg, ))
                                                                ), ), \
                        legs=sub.SubtractionLegSet((oth_leg, )))
                        if overall_parents[0] != None:

                            mapping_dict = bidict.bidict({overall_parents[0]: frozenset([singular_structure.legs[0].n, ia])})
#                            mapping_info_list.append(\
#                                ((ia,ib), [{'mapping_class': mapping_information['mapping'],
#                                            'mapping_singular_structure': mapping_structure,
#                                            'mapping_momenta_dict': mapping_dict}])\
#                                )

                        elif ia != ib:
#                        elif ia > ib:
                            mapping_dict = bidict.bidict({-1: frozenset([singular_structure.legs[0].n, sec_leg.n])})
#                            mapping_info_list.append(\
#                                ((ia,ib), [{'mapping_class': mapping_information['mapping'],
#                                            'mapping_singular_structure': mapping_structure,
#                                            'mapping_momenta_dict': mapping_dict}])\
#                                )
                        else:
                            continue
                    else:
                        continue

                    mapping_info_list.append(\
                        ((ia,ib), [{'mapping_class': mapping_information['mapping'],
                                    'mapping_singular_structure': mapping_structure,
                                    'mapping_momenta_dict': mapping_dict}])\
                        )
        #logger1.info('From DipoleCurrent : ' + str(mapping_info_list))
        return mapping_info_list


    def get_colored_legs(self, reduced_process, global_variables):
        """ returns a dictionary {leg_number : leg}
         with only colored legs in it"""

        c_legs_numbers = []
        #logger1.info('c_legs_numbers_pre_start : ' + str(c_legs_numbers))
        c_legs_numbers = self.get_colored_leg_numbers(reduced_process).keys()
        #logger1.info('c_legs_numbers_start : ' + str(c_legs_numbers))
        overall_parents = global_variables['overall_parents']
        #logger1.info('overall_parents : ' + str(overall_parents[0]))
        #logger1.info('c_legs_numbers_start_post : ' + str(c_legs_numbers))

        # necessary step since evaluating SC extra we have [a b c parent] as legs_number and not the initial-state child
        tmp_c_legs_numbers = c_legs_numbers
        for i in range(0,len(c_legs_numbers)):
            for j in range(0,len(c_legs_numbers)):
                if c_legs_numbers[i] == 2 and c_legs_numbers[j] == overall_parents[0]:
                    tmp_c_legs_numbers[j] = 1
                elif c_legs_numbers[i] == 1 and c_legs_numbers[j] == overall_parents[0]:
                    tmp_c_legs_numbers[j] = 2
                elif c_legs_numbers[i] == overall_parents[0] and c_legs_numbers[j] == 2:
                    tmp_c_legs_numbers[i] = 1
                elif c_legs_numbers[i] == overall_parents[0] and c_legs_numbers[j] == 1:
                    tmp_c_legs_numbers[i] = 2
                else:
                    continue

#        for i in range(0,len(c_legs_numbers)):
#            for j in range(0,len(c_legs_numbers)):
#                if c_legs_numbers[i] == 2 and c_legs_numbers[j] == overall_parents[0]:
#                    c_legs_numbers[j] = 1
#                elif c_legs_numbers[i] == 1 and c_legs_numbers[j] == overall_parents[0]:
#                    c_legs_numbers[j] = 2
#                elif c_legs_numbers[i] == overall_parents[0] and c_legs_numbers[j] == 2:
#                    c_legs_numbers[i] = 1
#                elif c_legs_numbers[i] == overall_parents[0] and c_legs_numbers[j] == 1:
#                    c_legs_numbers[i] = 2
#                else:
#                    continue
        #logger1.info('c_legs_numbers_end : ' + str(c_legs_numbers))
        #logger1.info('tmp_c_legs_numbers_end : ' + str(tmp_c_legs_numbers))

        c_legs = {}
        reduced_process_legs = []
        for l in reduced_process.get('legs'):
            reduced_process_legs.append(l)
        #logger1.info('reduced_process_legs : ' + str(reduced_process_legs))

        if reduced_process_legs[0]['number'] == overall_parents[0]:
            reduced_process_legs[0]['number'] = 1
        elif reduced_process_legs[1]['number'] == overall_parents[0]:
            reduced_process_legs[1]['number'] = 2
        #logger1.info('reduced_process_legs_post : ' + str(reduced_process_legs))


        for l in reduced_process.get('legs'):
            if l['number'] in tmp_c_legs_numbers:

                c_legs[l['number']] = sub.SubtractionLeg(l)

####
#        c_legs_numbers = self.get_colored_leg_numbers(reduced_process).keys()
#        if c_legs_numbers[0] == 6 and c_legs_numbers[1] == 1:
#            c_legs_numbers[0] = 2
#        elif c_legs_numbers[0] == 6 and c_legs_numbers[1] == 2:
#            c_legs_numbers[0] = 1
#        elif c_legs_numbers[0] == 1 and c_legs_numbers[1] == 6:
#            c_legs_numbers[1] = 2
#        elif c_legs_numbers[0] == 2 and c_legs_numbers[1] == 6:
#            c_legs_numbers[1] = 1
##        logger1.info('c_legs_numbers : ' + str(reduced_process_legs))

#        c_legs = {}
#        reduced_process_legs = []
#        for l in reduced_process.get('legs'):
#            reduced_process_legs.append(l)
##        logger1.info('c_legs_numbers : ' + str(reduced_process_legs[0]))


#        if reduced_process_legs[0]['number'] == 6:
#            reduced_process_legs[0]['number'] = 1
#        elif reduced_process_legs[1]['number'] == 6:
#            reduced_process_legs[1]['number'] = 2

#        for k in reduced_process_legs:
#            if k['number'] in c_legs_numbers:
#                #logger1.info('l number :' + str(l['number']))
##                if l['number'] == 6:
##                    l['number'] = 1

#                c_legs[k['number']] = sub.SubtractionLeg(k)

####

##        c_legs = {}
#        for l in reduced_process.get('legs'):
##            if l['number'] == 6:
##                l['number'] = 1
#            if l['number'] in c_legs_numbers:
#                #logger1.info('l number :' + str(l['number']))
##                if l['number'] == 6:
##                    l['number'] = 1

#                c_legs[l['number']] = sub.SubtractionLeg(l)


 
        #logger1.info('clegs :' + str(c_legs))

        return c_legs






