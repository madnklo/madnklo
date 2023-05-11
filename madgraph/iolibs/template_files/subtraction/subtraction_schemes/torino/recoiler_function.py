class MadEvent7Error(Exception):
    pass

"""Function for choosing the reference particle"""
" It returns the leg of the reference particle, i.d. tuple of (number, id, state) "

def get_recoiler(defining_process, excluded=()):

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

    # # for integrated currents
    # if len(sector_legs) == 0:
    #     if global_variables['overall_children'][0][0] > 2 and global_variables['overall_children'][0][1] > 2:   # final splitting
    #         if len(initial_recoilers) > 0:
    #             #print('caso F rec I : ' + str(sub.SubtractionLegSet([initial_recoilers[0]])))
    #             return sub.SubtractionLegSet([initial_recoilers[0]])
    #         else: 
    #             # sort the recoilers according to their id, and return the first one
    #             final_recoilers.sort(key = lambda l: l['id'])
    #             return sub.SubtractionLegSet([final_recoilers[0]])
    #     if global_variables['overall_children'][0][0] <= 2 or global_variables['overall_children'][0][1] <= 2:  # initial splitting
    #         if len(initial_recoilers) > 0:
    #             #print('caso I rec I : ' + str(sub.SubtractionLegSet([initial_recoilers[0]])))
    #             return sub.SubtractionLegSet([initial_recoilers[0]])
    #         else: 
    #             # sort the recoilers according to their id, and return the first one
    #             final_recoilers.sort(key = lambda l: l['id'])
    #             #print('caso I rec F : ' + str(sub.SubtractionLegSet([final_recoilers[0]])))
    #             return sub.SubtractionLegSet([final_recoilers[0]])


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


def get_virtual_recoiler(leg_PDGs_proc_prefix):

    partons_id = [1, -1, 2, -2, 3, -3, 4, -4, 5, -5, 6, -6, 21]

    virtual_recoiler = []

    for i in range(0,len(leg_PDGs_proc_prefix)):

        parent = i

        if leg_PDGs_proc_prefix[parent] not in partons_id:
            continue

        initial_recoilers = [
            leg for leg in range(0,len(leg_PDGs_proc_prefix)) if all([
                leg <= 1,
                leg_PDGs_proc_prefix[leg] in partons_id,
                leg != parent]
            )
        ]

        final_recoilers = [
            leg for leg in range(0,len(leg_PDGs_proc_prefix)) if all([
                leg > 1,
                leg_PDGs_proc_prefix[leg] in partons_id,
                leg != parent]
            )
        ]

        if len(initial_recoilers) > 0:
            ref = initial_recoilers[0]
        else:
            final_recoilers.sort(key = lambda l: leg_PDGs_proc_prefix[l])
            ref = final_recoilers[0]

        virtual_recoiler.append((parent+1,ref+1))

    return virtual_recoiler
        

# def get_collinear_mapped_labels(a,b,c,isec,jsec, R_leg_PDGs, B_leg_PDGs):

#     if a <= 2:
#         raise MadEvent7Error('The first particle must be in the final state! (%d, %d, %d)' % (a,b,c))

#     # Associate a,b,c to the collinear particles in the sector
#     if a == isec:
#         rm_leg = a
#         if b == jsec:
#             parent_leg = b
#         elif c == jsec:
#             parent_leg = b
#     elif a == jsec:
#         rm_leg = a
#         if b == isec:
#             parent_leg = b
#         elif c == isec:
#             parent_leg = b
#     elif (rm_leg == 0 or parent_leg == 0 or rm_leg == parent_leg):
#         raise MadEvent7Error('Mapping tuple (%d%d%d) does not include particles defining the singular sector (%d%d)!'
#                                 % (a,b,c,isec,jsec))

#     mapped_flavours = R_leg_PDGs

#     if jsec > 2:
#         # q > q + g
#         if R_leg_PDGs[rm_leg] != 21 and R_leg_PDGs[parent_leg] == 21:
#             mapped_flavours[parent_leg] = R_leg_PDGs(rm_leg)
#         # q > g + q
#         elif R_leg_PDGs[rm_leg] == 21 and R_leg_PDGs[parent_leg] != 21:
#             mapped_flavours[parent_leg] = R_leg_PDGs(parent_leg)
#         # g > q + qx
#         elif R_leg_PDGs[rm_leg] == (- R_leg_PDGs[parent_leg]):
#             mapped_flavours[parent_leg] = 21
#         # g > g + g
#         elif R_leg_PDGs[rm_leg] == 21 and R_leg_PDGs[parent_leg] == 21:
#             mapped_flavours[parent_leg] = 21

#         mapped_flavours[rm_leg] = 0

#         for i in range(1, len(R_leg_PDGs) - 1)
#             for j in range(1, len(R_leg_PDGs))
#                 if B_leg_PDGs[i] == mapped_flavours[j]:
#                     if mapp