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
        

def get_collinear_mapped_labels(a,b,c,isec,jsec, R_leg_PDGs, B_leg_PDGs):
    """Labels and flavours of the n-body (Born) kinematics obtained from the
    (n+1)-body one by the mapping (abc) of sector (isec,jsec): isec (always
    a final-state parton) is removed and jsec becomes the parent [ij].
    The mapping is (ijr) or, for jsec initial with a final-state recoiler,
    (irj) (2209.09123, eq. (eq:CijRdef)).
    Returns the flavours of the (n+1) partons after the clustering (0 for
    isec, the flavour of [ij] for jsec), the Born label of each (n+1)
    parton (that of jsec for isec) and the parent leg jsec. Initial-state
    partons are matched to initial-state Born slots, final-state ones to
    final-state slots, by flavour and in order (as in get_mapped_labels,
    imap.f)."""

    if a != isec or isec <= 2 or jsec not in (b, c):
        raise MadEvent7Error('Mapping (%d,%d,%d) does not match sector (%d,%d)!'
                             % (a, b, c, isec, jsec))

    rm_leg = isec
    parent_leg = jsec
    f_rm = R_leg_PDGs[rm_leg - 1]
    f_parent = R_leg_PDGs[parent_leg - 1]

    # flavour of the parent [ij]
    if jsec > 2:
        # final-state splitting [ij] -> i + j
        if f_rm == 21:
            f_ij = f_parent             # q > g + q, g > g + g
        elif f_parent == 21:
            f_ij = f_rm                 # q > q + g
        elif f_rm == -f_parent:
            f_ij = 21                   # g > q + qx
        else:
            raise MadEvent7Error('Inconsistent final-state pair (%d,%d)' % (f_rm, f_parent))
    else:
        # initial-state splitting j -> i + [ij], with i final
        if f_rm == 21:
            f_ij = f_parent             # q > g + q, g > g + g
        elif f_parent == 21:
            f_ij = -f_rm                # g > q + qx (qx enters the Born)
        elif f_rm == f_parent:
            f_ij = 21                   # q > q + g (g enters the Born)
        else:
            raise MadEvent7Error('Inconsistent initial-state pair (%d,%d)' % (f_rm, f_parent))

    mapped_flavours = list(R_leg_PDGs)
    mapped_flavours[parent_leg - 1] = f_ij
    mapped_flavours[rm_leg - 1] = 0

    mapped_labels = [0] * len(R_leg_PDGs)
    assigned = [False] * len(B_leg_PDGs)
    for j in range(len(R_leg_PDGs)):
        if j == rm_leg - 1:
            continue
        for i in range(len(B_leg_PDGs)):
            if assigned[i] or (i <= 1) != (j <= 1):
                continue
            if B_leg_PDGs[i] == mapped_flavours[j]:
                assigned[i] = True
                mapped_labels[j] = i + 1
                break
    if not all(assigned):
        raise MadEvent7Error('Labels of the real (%s) and Born (%s) do not match!'
                             % (R_leg_PDGs, B_leg_PDGs))
    mapped_labels[rm_leg - 1] = mapped_labels[parent_leg - 1]

    return mapped_flavours, mapped_labels, parent_leg
