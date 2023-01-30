"""Function for choosing the reference particle"""
" It returns the leg of the reference particle, i.d. tuple of (number, id, state) "

def get_recoiler(defining_process, excluded=()):

    #model = reduced_process.get('model')
    try:
        sector_legs = excluded
    except AttributeError:
        # this is in the case of the integrated CT's
        sector_legs = ()

    partons_id = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21]

    initial_recoilers = [
        leg for leg in defining_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['state'] == leg.INITIAL,
            leg['number'] not in sector_legs])
                ]

    final_recoilers = [
        leg for leg in defining_process.get('legs') if all([
            leg['id'] in partons_id,
            leg['state'] == leg.FINAL,
            leg['number'] not in sector_legs])
                ]

#    # for integrated currents
#    if len(sector_legs) == 0:
#        if global_variables['overall_children'][0][0] > 2 and global_variables['overall_children'][0][1] > 2:   # final splitting
#            if len(initial_recoilers) > 0:
#                #print('caso F rec I : ' + str(sub.SubtractionLegSet([initial_recoilers[0]])))
#                return sub.SubtractionLegSet([initial_recoilers[0]])
#            else: 
#                # sort the recoilers according to their id, and return the first one
#                final_recoilers.sort(key = lambda l: l['id'])
#                return sub.SubtractionLegSet([final_recoilers[0]])
#        if global_variables['overall_children'][0][0] <= 2 or global_variables['overall_children'][0][1] <= 2:  # initial splitting
#            if len(initial_recoilers) > 0:
#                #print('caso I rec I : ' + str(sub.SubtractionLegSet([initial_recoilers[0]])))
#                return sub.SubtractionLegSet([initial_recoilers[0]])
#            else: 
#                # sort the recoilers according to their id, and return the first one
#                final_recoilers.sort(key = lambda l: l['id'])
#                #print('caso I rec F : ' + str(sub.SubtractionLegSet([final_recoilers[0]])))
#                return sub.SubtractionLegSet([final_recoilers[0]])


    # for local currents
    if sector_legs[0] > 2 and sector_legs[1] > 2:   # final splitting
        if len(initial_recoilers) > 0:
            return initial_recoilers[0]
        else:
            # sort the recoilers according to their id, and return the first one
            final_recoilers.sort(key = lambda l: l['id'])
            return final_recoilers[0]
    elif sector_legs[0] <= 2 or sector_legs[1] <= 2:    # initial splitting
        if len(initial_recoilers) > 0:
            return initial_recoilers[0]
        else:
            # sort the recoilers according to their id, and return the first one
            final_recoilers.sort(key = lambda l: l['id'])
            return final_recoilers[0]

