###########################################################
#
# Torino subtraction scheme -- Dipole Current class
#
###########################################################

import commons.general_current as general_current
import madgraph.core.subtraction as sub


class DipoleCurrent(general_current.GeneralCurrent):
    """A new daughter class which makes it possible to
    specify different recoilers for different color dipoles
    """

    def get_recoiler_combinations(self, reduced_process, global_variables, *args, **opts):
        """ This function returns a list of tuples of the following format:
             ( reduced_kinematic_identifier, (
                                        recoilers_subtraction_leg_set_for_step_1,
                                        recoilers_subtraction_leg_set_for_step_2,
                                        ...
                                        )
            )
            By default this returns ( None, []) and the list of recoiler legs from the entries 'reduced_recoilers'
            and 'additional_recoilers' of the mapping rules.
            Note that the leg numbers of these set of recoilers specified should be already mapped, i.e correspond
            to the leg numbers of the actual PS point and process supplied at the time this current is called.

            For Catani-Seymour dipole subtraction for example, you would use a different mapping for each dipole, so one
            could return a list of entries like the following for the dipole with legs (1,2):
            (   (1,2),   [  SubtractionLegSet(SubtractionLeg(number=2)), ]  )
        """

        recoiler_list = []

        for ia, lega in self.get_colored_legs(reduced_process).items():
            for ib, legb in self.get_colored_legs(reduced_process).items():
                if ia == ib:
                    continue
                recoiler_list.append(\
                    ((ia,ib), [sub.SubtractionLegSet((lega, legb)),])\
                    )
        return recoiler_list


    def get_colored_legs(self, reduced_process):
        """ returns a dictionary {leg_number : leg}
         with only colored legs in it"""

        c_legs_numbers = self.get_colored_leg_numbers(reduced_process).keys()
        c_legs = {}
        for l in reduced_process.get('legs'):
            if l['number'] in c_legs_numbers:
                c_legs[l['number']] = sub.SubtractionLeg(l)

        return c_legs






