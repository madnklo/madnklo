################################################################################
#
# Copyright (c) 2009 The MadGraph5_aMC@NLO Development team and Contributors
#
# This file is a part of the MadGraph5_aMC@NLO project, an application which 
# automatically generates Feynman diagrams and matrix elements for arbitrary
# high-energy processes in the Standard Model and beyond.
#
# It is subject to the MadGraph5_aMC@NLO license which should accompany this 
# distribution.
#
# For more information, visit madgraph.phys.ucl.ac.be and amcatnlo.web.cern.ch
#
################################################################################

"""Classes, methods and functions required to write QCD color information 
for a diagram and build a color basis, and to square a QCD color string for
squared diagrams and interference terms."""

import copy
import fractions
import operator
import re
import array

import madgraph.core.color_algebra as color_algebra
import madgraph.core.diagram_generation as diagram_generation
import madgraph.core.base_objects as base_objects
from madgraph import MadGraph5Error
import madgraph.various.misc as misc

#===============================================================================
# ColorBasis
#===============================================================================
class ColorBasis(dict):
    """The ColorBasis object is a dictionary created from an amplitude. Keys
    are the different color structures present in the amplitude. Values have
    the format (diag,(index c1, index c2,...), coeff, is_imaginary, Nc_power) 
    where diag is the diagram index, (index c1, index c2,...) the list of 
    indices corresponding to the chose color parts for each vertex in the 
    diagram, coeff the corresponding coefficient (a fraction), is_imaginary
    if this contribution is real or complex, and Nc_power the Nc power."""

    # Dictionary to save simplifications already done in a canonical form
    _canonical_dict = {}

    # Dictionary store the raw colorize information
    _list_color_dict = []


    class ColorBasisError(Exception):
        """Exception raised if an error occurs in the definition
        or the execution of a color basis object."""
        pass

    def colorize(self, diagram, model):
        """Takes a diagram and a model and outputs a dictionary with keys being
        color coefficient index tuples and values a color string (before 
        simplification)."""

        # The smallest value used to create new summed indices
        min_index = -1000
        # The dictionary to be output
        res_dict = {}
        # The dictionary for book keeping of replaced indices
        repl_dict = {}

        for i, vertex in enumerate(diagram.get('vertices')):
            min_index, res_dict = self.add_vertex(vertex, diagram, model,
                            repl_dict, res_dict, min_index)

        # if the process has no QCD particles
        # Return a list filled with ColorOne if all entries are empty ColorString()
        empty_colorstring = color_algebra.ColorString()
        if all(cs == empty_colorstring for cs in res_dict.values()):
            res_dict = dict((key, color_algebra.ColorString(
                               [color_algebra.ColorOne()])) for key in res_dict)
                    
        return res_dict

    

    def add_vertex(self, vertex, diagram, model,
                   repl_dict, res_dict, min_index, id0_rep=[]):
        """Update repl_dict, res_dict and min_index for normal vertices.
        Returns the min_index reached and the result dictionary in a tuple.
        If the id0_rep list is not None, perform the requested replacement on the
        last leg number before going further."""

        # Create a list of (color,leg number) pairs for the vertex, where color
        # can be negative for anti particles

        color_num_pairs = []
        pdg_codes = []
                
        for index, leg in enumerate(vertex.get('legs')):
            curr_num = leg.get('number')
            curr_part = model.get('particle_dict')[leg.get('id')]
            curr_color = curr_part.get_color()
            curr_pdg = curr_part.get_pdg_code()

            # If this is the next-to-last vertex and the last vertex is
            # the special identity id=0, start by applying the replacement rule
            # on the last vertex.
            if index == len(vertex.get('legs')) - 1 and \
                    curr_num in id0_rep:
                    curr_num = id0_rep[id0_rep.index(curr_num) - 1]

            # If this is the last leg and not the last vertex 
            # flip color. If it is not the last, AND not the next-to-last
            # before an id=0 vertex, replace last index by a new summed index.
            if index == len(vertex.get('legs')) - 1 and \
                vertex != diagram.get('vertices')[-1]:
                curr_color = curr_part.get_anti_color()
                curr_pdg = curr_part.get_anti_pdg_code()
                if not id0_rep:
                    if not ( diagram.get('vertices')[-1].get('id')==-1 and \
                    vertex == diagram.get('vertices')[-2]):
                        repl_dict[curr_num] = min_index
                        min_index = min_index - 1
                    else:                  
                        repl_dict[curr_num] = \
                          max(l.get('number') for l in \
                                        diagram.get('vertices')[-1].get('legs'))

            # Take into account previous replacements
            try:
                curr_num = repl_dict[curr_num]
            except KeyError:
                pass

            color_num_pairs.append((curr_color, curr_num))
            pdg_codes.append(curr_pdg)

        if vertex != diagram.get('vertices')[-1]:
            # Put the resulting wavefunction first, to make
            # wavefunction call more natural
            last_color_num = color_num_pairs.pop(-1)
            color_num_pairs.insert(0, last_color_num)
            last_pdg = pdg_codes.pop(-1)
            pdg_codes.insert(0, last_pdg)

        # Order the legs according to the interaction particles
        if vertex.get('id')!=-1:
            interaction_pdgs = [p.get_pdg_code() for p in \
                                model.get_interaction(vertex.get('id')).\
                                get('particles')]
        else:
            interaction_pdgs = [l.get('id') for l in vertex.get('legs')]

        sorted_color_num_pairs = []
        #print "interactions_pdg=",interaction_pdgs
        #print "pdg_codes=",pdg_codes        
        for i, pdg in enumerate(interaction_pdgs):
            index = pdg_codes.index(pdg)
            pdg_codes.pop(index)
            sorted_color_num_pairs.append(color_num_pairs.pop(index))

        if color_num_pairs:
            raise base_objects.PhysicsObject.PhysicsObjectError

        color_num_pairs = sorted_color_num_pairs

        # Create a list of associated leg number following the same order
        list_numbers = [p[1] for p in color_num_pairs]

        # ... and the associated dictionary for replacement
        match_dict = dict(enumerate(list_numbers))

        if vertex['id'] == -1:
            return (min_index, res_dict)

        # Update the result dict using the current vertex ColorString object
        # If more than one, create different entries
        inter_color = model.get_interaction(vertex['id'])['color']
        inter_indices = [i for (i,j) in \
                        model.get_interaction(vertex['id'])['couplings'].keys()]
        
        # For colorless vertices, return a copy of res_dict
        # Where one 0 has been added to each color index chain key
        if not inter_color:
            new_dict = {}
            for k, v in res_dict.items():
                new_key = tuple(list(k) + [0])
                new_dict[new_key] = v
            # If there is no result until now, create an empty CS...
            if not new_dict:
                new_dict[(0,)] = color_algebra.ColorString()
            return (min_index, new_dict)

        new_res_dict = {}
        for i, col_str in \
                enumerate(inter_color):
            
            # Ignore color string if it doesn't correspond to any coupling
            if i not in inter_indices:
                continue
            
            # Build the new element
            assert type(col_str) == color_algebra.ColorString 
            mod_col_str = col_str.create_copy()

            # Replace summed (negative) internal indices
            list_neg = []
            for col_obj in mod_col_str:
                list_neg.extend([ind for ind in col_obj if ind < 0])
            internal_indices_dict = {}
            # This notation is to remove duplicates
            for index in list(set(list_neg)):
                internal_indices_dict[index] = min_index
                min_index = min_index - 1
            mod_col_str.replace_indices(internal_indices_dict)

            # Replace other (positive) indices using the match_dic
            mod_col_str.replace_indices(match_dict)

            # If we are considering the first vertex, simply create
            # new entries

            if not res_dict:
                new_res_dict[tuple([i])] = mod_col_str
            #... otherwise, loop over existing elements and multiply
            # the color strings
            else:
                for ind_chain, col_str_chain in res_dict.items():
                    new_col_str_chain = col_str_chain.create_copy()
                    new_col_str_chain.product(mod_col_str)
                    new_res_dict[tuple(list(ind_chain) + [i])] = \
                        new_col_str_chain

        return (min_index, new_res_dict)


    def update_color_basis(self, colorize_dict, index):
        """Update the current color basis by adding information from 
        the colorize dictionary (produced by the colorize routine)
        associated to diagram with index index. Keep track of simplification
        results for maximal optimization."""
        import madgraph.various.misc as misc
        # loop over possible color chains
        for col_chain, col_str in colorize_dict.items():
            # Create a canonical immutable representation of the the string
            canonical_rep, rep_dict = col_str.to_canonical()
            try:
                # If this representation has already been considered,
                # recycle the result.                               
                col_fact = self._canonical_dict[canonical_rep].create_copy()
            except KeyError:
                # If the representation is really new

                # Create and simplify a color factor for the considered chain
                col_fact = color_algebra.ColorFactor([col_str])
                col_fact = col_fact.full_simplify()

                # Here we need to force a specific order for the summed indices
                # in case we have K6 or K6bar Clebsch Gordan coefficients
                for colstr in col_fact: colstr.order_summation()

                # Save the result for further use
                canonical_col_fact = col_fact.create_copy()
                canonical_col_fact.replace_indices(rep_dict)
                # Remove overall coefficient
                for cs in canonical_col_fact:
                    cs.coeff = cs.coeff / col_str.coeff
                self._canonical_dict[canonical_rep] = canonical_col_fact
            else:
                # If this representation has already been considered,
                # adapt the result
                # Note that we have to replace back
                # the indices to match the initial convention. 
                col_fact.replace_indices(self._invert_dict(rep_dict))
                # Since the initial coeff of col_str is not taken into account
                # for matching, we have to multiply col_fact by it.
                for cs in col_fact:
                    cs.coeff = cs.coeff * col_str.coeff
                # Must simplify up to two times at NLO (since up to two traces
                # can appear with a loop) to put traces in a canonical ordering.
                # If it still causes issue, just do a full_simplify(), it would
                # not bring any heavy additional computational load.
                col_fact = col_fact.simplify().simplify()
                
                # Here we need to force a specific order for the summed indices
                # in case we have K6 or K6bar Clebsch Gordan coefficients
                for colstr in col_fact: colstr.order_summation()

            # loop over color strings in the resulting color factor
            for col_str in col_fact:
                immutable_col_str = col_str.to_immutable()
                # if the color structure is already present in the present basis
                # update it
                basis_entry = (index,
                                col_chain,
                                col_str.coeff,
                                col_str.is_imaginary,
                                col_str.Nc_power,
                                col_str.loop_Nc_power)
                try:
                    self[immutable_col_str].append(basis_entry)
                except KeyError:
                    self[immutable_col_str] = [basis_entry]

    def create_color_dict_list(self, amplitude):
        """Returns a list of colorize dict for all diagrams in amplitude. Also
        update the _list_color_dict object accordingly """

        list_color_dict = []

        for diagram in amplitude.get('diagrams'):
            colorize_dict = self.colorize(diagram,
                                          amplitude.get('process').get('model'))
            list_color_dict.append(colorize_dict)

        self._list_color_dict = list_color_dict

        return list_color_dict

    def build(self, amplitude=None):
        """Build the a color basis object using information contained in
        amplitude (otherwise use info from _list_color_dict). 
        Returns a list of color """

        if amplitude:
            self.create_color_dict_list(amplitude)
        for index, color_dict in enumerate(self._list_color_dict):
            self.update_color_basis(color_dict, index)

    def __init__(self, *args):
        """Initialize a new color basis object, either empty or filled (0
        or 1 arguments). If one arguments is given, it's interpreted as 
        an amplitude."""

        assert len(args) < 2, "Object ColorBasis must be initialized with 0 or 1 arguments"


        dict.__init__(self)

        # Dictionary to save simplifications already done in a canonical form
        self._canonical_dict = {}

        # Dictionary store the raw colorize information
        self._list_color_dict = []


        if args:
            assert isinstance(args[0], diagram_generation.Amplitude), \
                        "%s is not a valid Amplitude object" % str(args[0])
                        
            self.build(*args)

    def __str__(self):
        """Returns a nicely formatted string for display"""

        my_str = ""
        for k, v in self.items():
            for name, indices in k:
                my_str = my_str + name + str(indices)
            my_str = my_str + ': '
            for contrib in v:
                imag_str = ''
                if contrib[3]:
                    imag_str = 'I'
                my_str = my_str + '(diag:%i, chain:%s, coeff:%s%s, Nc:%i) ' % \
                                    (contrib[0], contrib[1], contrib[2],
                                     imag_str, contrib[4])
            my_str = my_str + '\n'
        return my_str

    def _invert_dict(self, mydict):
        """Helper method to invert dictionary dict"""

        return dict([v, k] for k, v in mydict.items())

    @staticmethod
    def get_color_flow_string(my_color_string, octet_indices):
        """Return the color_flow_string (i.e., composed only of T's with 2 
        indices) associated to my_color_string. Take a list of the external leg
        color octet state indices as an input. Returns only the leading N 
        contribution!"""
        # Create a new color factor to allow for simplification
        my_cf = color_algebra.ColorFactor([my_color_string])

        # Add one T per external octet
        for indices in octet_indices:
            if indices[0] == -6:
                # Add a K6 which contracts the antisextet index to a
                # pair of antitriplets
                my_cf[0].append(color_algebra.K6(indices[1],
                                                 indices[2],
                                                 indices[3]))
            if indices[0] == 6:
                # Add a K6Bar which contracts the sextet index to a
                # pair of triplets
                my_cf[0].append(color_algebra.K6Bar(indices[1],
                                                    indices[2],
                                                    indices[3]))
            if abs(indices[0]) == 8:
                # Add a T which contracts the octet to a
                # triplet-antitriplet pair
                my_cf[0].append(color_algebra.T(indices[1],
                                                indices[2],
                                                indices[3]))
        # Simplify the whole thing
        my_cf = my_cf.full_simplify()

        # If the result is empty, just return
        if not my_cf:
            return my_cf

        # Return the string with the highest N coefficient 
        # (leading N decomposition), and the value of this coeff
        max_coeff = max([cs.Nc_power for cs in my_cf])

        res_cs = [cs for cs in my_cf if cs.Nc_power == max_coeff]

        # If more than one string at leading N...
        if len(res_cs) > 1 and any([not cs.near_equivalent(res_cs[0]) \
                                    for cs in res_cs]):
            raise ColorBasis.ColorBasisError, \
             "More than one color string with leading N coeff: %s" % str(res_cs)

        res_cs = res_cs[0]

        # If the result string does not contain only T's with two indices
        # and Epsilon/EpsilonBar objects
        for col_obj in res_cs:
            if not isinstance(col_obj, color_algebra.T) and \
                   not col_obj.__class__.__name__.startswith('Epsilon'):
                raise ColorBasis.ColorBasisError, \
                  "Color flow decomposition %s contains non T/Epsilon elements" % \
                                                                    str(res_cs)
            if isinstance(col_obj, color_algebra.T) and len(col_obj) != 2:
                raise ColorBasis.ColorBasisError, \
                  "Color flow decomposition %s contains T's w/o 2 indices" % \
                                                                    str(res_cs)

        return res_cs

    def color_flow_decomposition(self, repr_dict, ninitial):
        """Returns the color flow decomposition of the current basis, i.e. a 
        list of dictionaries (one per color basis entry) with keys corresponding
        to external leg numbers and values tuples containing two color indices
        ( (0,0) for singlets, (X,0) for triplet, (0,X) for antitriplet and 
        (X,Y) for octets). Other color representations are not yet supported 
        here (an error is raised). Needs a dictionary with keys being external
        leg numbers, and value the corresponding color representation."""

        # Offsets used to introduce fake quark indices for gluons
        offset1 = 1000
        offset2 = 2000
        offset3 = 3000

        res = []

        for col_basis_entry in sorted(self.keys()):

            res_dict = {}
            fake_repl = []

            # Rebuild a color string from a CB entry
            col_str = color_algebra.ColorString()
            col_str.from_immutable(col_basis_entry)
            for (leg_num, leg_repr) in repr_dict.items():
                # By default, assign a (0,0) color flow
                res_dict[leg_num] = [0, 0]

                # Raise an error if external legs contain non supported repr
                if abs(leg_repr) not in [1, 3, 6, 8]:
                    raise ColorBasis.ColorBasisError, \
        "Particle ID=%i has an unsupported color representation" % leg_repr

                # Build the fake indices replacements for octets
                if abs(leg_repr) == 8:
                    fake_repl.append((leg_repr, leg_num,
                                      offset1 + leg_num,
                                      offset2 + leg_num))
                # Build the fake indices for sextets
                elif leg_repr in [-6, 6]:
                    fake_repl.append((leg_repr, leg_num,
                                      offset1 + leg_num,
                                      offset3 + leg_num))

            # Get the actual color flow
            col_str_flow = self.get_color_flow_string(col_str, fake_repl)

            # Offset for color flow
            offset = 500

            for col_obj in col_str_flow:
                if isinstance(col_obj, color_algebra.T):
                    # For T, all color indices should be the same
                    offset = offset + 1
                for i, index in enumerate(col_obj):
                    if isinstance(col_obj, color_algebra.Epsilon):
                        # Epsilon contracts with antitriplets,
                        i = 0
                        # ...and requires all different color indices
                        offset = offset+1
                    elif isinstance(col_obj, color_algebra.EpsilonBar):
                        # EpsilonBar contracts with antitriplets
                        i = 1
                        # ...and requires all different color indices
                        offset = offset+1
                    if index < offset1:
                        res_dict[index][i] = offset
                    elif index > offset1 and index < offset2:
                        res_dict[index - offset1][i] = offset
                    elif index > offset2 and index < offset3:
                        res_dict[index - offset2][i] = offset
                    elif index > offset3:
                        # For color sextets, use negative triplet
                        # number to reperesent antitriplet and vice
                        # versa, allowing for two triplet or two
                        # antitriplet numbers representing the color
                        # sextet.
                        res_dict[index - offset3][1-i] = -offset

            # Reverse ordering for initial state to stick to the (weird)
            # les houches convention

            for key in res_dict.keys():
                if key <= ninitial:
                    res_dict[key].reverse()

            res.append(res_dict)

        return res




#===============================================================================
# ColorMatrix
#===============================================================================
class ColorMatrix(dict):
    """A color matrix, meaning a dictionary with pairs (i,j) as keys where i
    and j refer to elements of color basis objects. Values are Color Factor
    objects. Also contains two additional dictionaries, one with the fixed Nc
    representation of the matrix, and the other one with the "inverted" matrix,
    i.e. a dictionary where keys are values of the color matrix."""

    # See doc of 'generate_all_color_connections' for details on the meaning of
    # these general parameters
    Q1_DUMMY_STARTING_INDEX               = -100000  
    Q2_DUMMY_STARTING_INDEX               = -200000
    SINGLE_PRIME_POSITIVE_INDICES_OFFSET  = 100 
    DOUBLE_PRIME_POSITIVE_INDICES_OFFSET  = 200
    EMITTED_NEGATIVE_INDICES_OFFSET       = 300

    def __init__(self, col_basis, col_basis2=None,
                 Nc=3, Nc_power_min=None, Nc_power_max=None, automatic_build=True):
        """Initialize a color matrix with one or two color basis objects. If
        only one color basis is given, the other one is assumed to be equal.
        As options, any value of Nc and minimal/maximal power of Nc can also be 
        provided. Note that the min/max power constraint is applied
        only at the end, so that it does NOT speed up the calculation."""

        self.col_matrix_fixed_Nc = {}
        self.inverted_col_matrix = {}
        
        self._col_basis1 = col_basis
        if col_basis2:
            self._col_basis2 = col_basis2
            self.is_symmetric = False
            if automatic_build:
                self.build_matrix(Nc, Nc_power_min, Nc_power_max)
        else:
            self._col_basis2 = col_basis
            self.is_symmetric = True
            if automatic_build:
                # If the two color basis are equal, assumes the color matrix is symmetric
                self.build_matrix(Nc, Nc_power_min, Nc_power_max, is_symmetric=True)

    def get_splitting_color_operator(self, incoming_index, emitting_repr, outgoing_base_index, 
                                             outgoing_emitted_index, qqbar=False, Q='Q1'):
        """ Returns the color operator associated with the splitting of emitting_index (in the
        emitting_repr representation (in [8,3,-3]) into the emitted_index.
        The flag qqbar indicates if this is a g > g g or g > q q~ splitting.
        The outgoing_emitted_index is the new negative emitted index; that is:
        
        The convention for T in the color module is T^{octet}_{incoming triplet, outgoing triplet}
        
          q > q g    :  T^{outgoing_emitted_index}_{incoming_index,outgoing_base_index}
          q~ > q~ g  :  -T^{outgoing_emitted_index}_{outgoing_base_index, incoming_index}
          g > g g    :  f^{incoming_index, outgoing_emitted_index, outgoing_base_index}
          g > q q~   :  T^{incoming_index}_{outgoing_base_index,outgoing_emitted_index}
          
        Finally, when this splitting is to be added to the Q2 colour structure, the convention
        is to *never* apply complex conjugation except for the g > q q~ splitting, where
        it is necessary because the emitted particle (and antiquark) is not its own self-antiparticle.
        In this case, the splitting operator added is the complex-conjugate of the one added
        for Q1:

          g > q q~ (for Q2) :  -T^{incoming_index}_{outgoing_emitted_index,outgoing_base_index}        
        
        """
        
        if emitting_repr not in [3,-3,8]:
            raise MadGraph5Error("Color connected matrix elements between particles color-charged"+
                                 " in the %d representation are not implemented."%abs(emitting_repr))
        if qqbar and emitting_repr != 8:
            raise MadGraph5Error("Only a color octet can split into a q-qbar pair.")

        # q > q g
        if emitting_repr == 3:
            return color_algebra.ColorString([color_algebra.T(outgoing_emitted_index, incoming_index, outgoing_base_index)])
        # q~ > q~ g
        elif emitting_repr == -3:
            # The anti quark emissiong operator color operator is (T^c_{ab})* = -T^c_{ba}
            return color_algebra.ColorString([color_algebra.T(outgoing_emitted_index, outgoing_base_index, incoming_index)],
                                                                                    coeff=fractions.Fraction(-1, 1))
        # For gluon self-interactions, we chose the second index 'b' of f^{abc} to be the one carrying 
        # "emitted" gluon's color.
        # The dual representation of the color charge takes an imaginary factor 'i', this is also what insures that
        # the color matrix remains real for the correlator of the type f * T
        elif (not qqbar) and emitting_repr == 8:
            return color_algebra.ColorString([color_algebra.f(incoming_index, outgoing_emitted_index, outgoing_base_index)],
                                                                                                 is_imaginary=True)
        # The outgoing_emitted_index is always chosen to be an outgoing anti-quark while the
        # outgoing_base_index is always a quark. 
        elif qqbar and emitting_repr == 8:
            if Q=='Q1' or True:
                return color_algebra.ColorString([color_algebra.T(incoming_index, 
                                            outgoing_emitted_index, outgoing_base_index)])
            else:
                return color_algebra.ColorString([color_algebra.T(incoming_index,
                    outgoing_base_index, outgoing_emitted_index)],coeff=fractions.Fraction(-1, 1))

    def add_splitting_to_connection(self, connection, emitting_index, emitting_repr, 
                                                              emitted_index, qqbar=False ):
        """ Add a splitting to the connection passed in argument. The splitting is specified
        by its emitting index and its representation as well as the emitted index.
        A special flag denotes if a g > g g or g > q q~ must be considered."""
        
        # First remove the emitted_index from the list of those that must still be emitted
        connection['emitted_numbers_pool'].remove(emitted_index)

        # Now add the correct color operator to the current color strings representing 
        # the connection so far. 
        # We start by determining the first outgoing index for Q1 and Q2
        if emitting_index > 0:
            first_outgoing_index_Q1  = self.SINGLE_PRIME_POSITIVE_INDICES_OFFSET+emitting_index
        else:
            first_outgoing_index_Q1 = self.EMITTED_NEGATIVE_INDICES_OFFSET+abs(emitting_index)
        if emitting_index > 0:
            first_outgoing_index_Q2  = self.DOUBLE_PRIME_POSITIVE_INDICES_OFFSET+emitting_index
        else:
            first_outgoing_index_Q2 = self.EMITTED_NEGATIVE_INDICES_OFFSET+abs(emitting_index)
        
        # The second outgoing index is always identical for Q1 and Q2 as it must be negative
        assert(emitted_index<0)
        second_outgoing_index = self.EMITTED_NEGATIVE_INDICES_OFFSET+abs(emitted_index)

        # We must now determine the emitting color index and replace it with a dummy summation
        # if needed. We start with Q1
        if emitting_index in connection['replace_indices_Q1']:
            # This operator will directly link to the matrix element, no dummy indext must
            # be introduced then
            del connection['replace_indices_Q1'][emitting_index]
            incoming_index_Q1 = emitting_index
        else:
            # We must introduce a dummy index
            connection['last_dummy_index_Q1'] -= 1
            incoming_index_Q1 = connection['last_dummy_index_Q1']
            # We must replace the exiting first_outgoing_index_Q1 in the colour structure
            # by this new dummy index
            connection['color_string_Q1'].replace_indices({first_outgoing_index_Q1: incoming_index_Q1})

        # And similarly for Q2
        if self.SINGLE_PRIME_POSITIVE_INDICES_OFFSET+emitting_index in connection['replace_indices_Q2']:
            # This operator will directly link to the matrix element, no dummy indext must
            # be introduced then
            del connection['replace_indices_Q2'][self.SINGLE_PRIME_POSITIVE_INDICES_OFFSET+emitting_index]
            incoming_index_Q2 = self.SINGLE_PRIME_POSITIVE_INDICES_OFFSET+emitting_index
        else:
            # We must introduce a dummy index
            connection['last_dummy_index_Q2'] -= 1
            incoming_index_Q2 = connection['last_dummy_index_Q2']
            
            # We must replace the exiting first_outgoing_index_Q1 in the colour structure
            # by this new dummy index
            connection['color_string_Q2'].replace_indices({first_outgoing_index_Q2: incoming_index_Q2})
            
        # Now that all the indices of the color operator to add have been figured out,
        # we can add it to the corresponding color strings
        connection['color_string_Q1'].product( self.get_splitting_color_operator(
            incoming_index_Q1, emitting_repr, first_outgoing_index_Q1, second_outgoing_index, qqbar=qqbar, Q='Q1') )
        connection['color_string_Q2'].product( self.get_splitting_color_operator(
            incoming_index_Q2, emitting_repr, first_outgoing_index_Q2, second_outgoing_index, qqbar=qqbar, Q='Q2') )
    
        # We must add the newly generated index to the list of available ones. 
        if not qqbar:
            # It is always a gluon then
            connection['open_indices'].append((emitted_index,8))
            # Add the splitting to the tuple representation
            connection['tuple_representation'].append( 
                                             (emitting_index, emitted_index, emitting_index) )
        else:
            # We emit an anti-quark in this case
            connection['open_indices'].append((emitted_index,-3))
            # We must remove the existing open index which was a gluon and replace it with
            # a quark one
            connection['open_indices'].remove((emitting_index,8))
            connection['open_indices'].append((emitting_index,3))
            # As you can see the convention for such splitting is to be denoted 
            # (-1,-1,-2) while g > g g would be labelled as (-1,-2,-1).
            connection['tuple_representation'].append( 
                                          (emitting_index, emitting_index, emitted_index) )
            

    def generate_all_color_connections(self, process_legs, model, order='NLO'):
        """Returns a dictionary whose keys is the "identifier" of the color connection.
        This identifier is formatted as the following 2-tuple structure:
            ( first_connection_structure, second_connection_structure )
        where each connection structure is a tuple of tuples formatted as follows:
            ( ColorOperator1, ColorOperator2, ColorOperator3, etc.. )
        where each color operator (Q) is a tuple with the following three entries:
            ( incoming_saturated_index (i), new_open_index_A (j), new_open_index_B (k) )
        which identifies the color indices "generated" by the color operator which is 
        connected to the previously open index labeled 'incoming_saturated_index'.
        
        Depending of the color representation of the 'incoming_saturated_index', one has:
          A) color-triplet      : Q -> T^j_{ik}
          B) color-anti-triplet : Q -> -T^j_{ik}
          C) color-octet        : Q -> i f^ijk or
          D) color-octet        : Q -> T^i_{jk} if i == j and i<0
        The last condition i<0 enforces that the g > q q~ splitting occurs from an already
        soft gluon.
        
        We now illustrate the notation above with the process e+(1) e-(2) > d(3) d~(4) g(5):
        At NLO the connections are simply the following three:
            ((3,-1,3),) 
            ((4,-1,4),) 
            ((5,-1,5),)
        Leading to 9 possible color correlators. At NNLO, the possible soft currents are:
        (Notice that using our conventions, the list of color operators is unordered and
        can therefore be considered as a set. For convenience however, we will considered
        its canonical representation to be an ordered list instead. The examples below
        do not reflect this ordering and instead you can imagine all of them being wrapped
        by a "sorted(*, key= lambda el: el[0], reverse=True)" statement. Notice that because
        of the dummy indices, the resulting basis might still not be minimal).
        
          > The "simple disjoint dipoles"
              ( (3,-1,3), (4,-2,4) )
              ( (3,-1,3), (5,-2,5) )
              ( (4,-1,4), (5,-2,5) )
            
          > The "successive gluon emissions from a single quark leg"
              ( (3,-1,3), (3,-2,3) )
              ( (4,-1,4), (4,-2,4) )
            
          > The "successive gluon emissions from a single gluon leg"
              ( (5,-1,5), (5,-2,5) )
              ( (5,-1,5), (-1,-2,-1) )
          
          > The "gluon splitting into a quark pair" 
              ( (3,-1,3), (-1,-1,-2) )
              ( (4,-1,4), (-1,-1,-2) )
              ( (5,-1,5), (-1,-1,-2) )
            Notice how the (5,-1,-2) operator is not allowed and how (-1,-1,-2) indicates
            the splitting into two quarks while (-1,-2,-1) before indicated the splitting
            into two gluons.

          >  Emission from one quark and the gluon
            ( (3,-1,3), (5,-2,5) )
            ( (4,-1,4), (5,-2,5) )
            
        Now, we add that each of the color current structure above will exist in two copies
        with the label -1 and -2 interchanged *when not repeated/summed over* and *when not part of
        the same tuple*. So '( (5,-1,5), (-1,-2,-1) )' does not have any symmetric counterpart.
        An analoguous symmetrisation involving k gluon indices is performed in an N^kLO computation.
        

        Finally, and more pragmatically, the computation is performed as follows:
        
        < M(c_i) | Q1(c_i -> c_i') Q2(c_i' -> c_i'') | M(c_i'') >
        
        where c_i are sets of color indices and the "primes" ' indicate the output images
        of the action of color operators Q. 
        The dummy indices part of < M(c_i) | and | M(c_i'') > will be forced to be non-overlapping,
        while those for Q1(c_i -> c_i') (resp. Q2) will be selected decrementally from
        Q1_DUMMY_STARTING_INDEX (resp. Q1_DUMMY_STARTING_INDEX).
        The set of indices c_i are positive and span all leg numbers.
        The set of color indices c_i' and c_i'' is constructed as follows:
        
            Negative indices c = (-1, -2, ...) which used to indicate gluon or quark emission
            are replaced by EMITTED_NEGATIVE_INDICES_OFFSET + |c| and positive ones
            by SINGLE_PRIME_POSITIVE_INDICES_OFFSET + c.
            
            Similarly for the c_i'' indices, using the variables DOUBLE_PRIME_POSITIVE_INDICES_OFFSET
            for positive ones and EMITTED_NEGATIVE_INDICES_OFFSET again for negative ones.
        
        We illustrate below the actual implementation of some of the correlators shown before.
        We will use here:
        
                Q1_DUMMY_STARTING_INDEX               = -100000  
                Q2_DUMMY_STARTING_INDEX               = -200000
                SINGLE_PRIME_POSITIVE_INDICES_OFFSET  = 100 
                DOUBLE_PRIME_POSITIVE_INDICES_OFFSET  = 200
                EMITTED_NEGATIVE_INDICES_OFFSET       = 300
                
    
        And the amplitude e+(1) e-(2) > d(3) d~(4) g(5) is denoted < 3 4 5 |, then the
        color connected ME for the connection ( (3,-1,3), (4,-2,4) ) squared against itself 
        will read:
        
        < 3 4 5 | 
                   T^{301}_{3,103} T^{302}_{4,104} \delta_{5,105}  
                   T^{301}_{103,203} T^{302}_{104,204} \delta_{105,205}
        | 203 204 205 >
        
        The connection ( (3,-1,3), (4,-2,4) ) squared against ( (3,-2,3), (5,-1,5) ) reads:

        < 3 4 5 | T^{301}_{3,103} T^{302}_{4,104} \delta_{5,105} 
                   \delta_{104,204} T^{302}_{103,203} i f^{105,301,205} | 203 204 205 >   
        
        For ( (3,-1,3), (4,-2,4) ) squared against ( (5,-1,5), (-1,-2,-1) ) we have:
        
        < 3 4 5 | T^{301}_{3,103} T^{302}_{4,104} \delta_{3,103} \delta_{4,104} 
                  \delta_{103,203} \delta_{104,204} i f^{105,-200001,205} i f^{-200001,302,301}
        | 203 204 205 >   
        
        and at last, for ( (5,-1,5), (-1,-1,-2) ) squared against itself:
        
        < 3 4 5 | \delta_{3,103} \delta_{4,104} i f^{5,-100001,105} T^{-100001}_{301,302}
                  \delta_{103,203} \delta_{104,204} i f^{105,-200001,205} T^{-200001}_{301,302}
        | 203 204 205 >   
        
        """
        
        max_order = order.count('N')
        
        
        # First collect all available external color indices and store them in a map with
        # which acts as the "seed" LO color correlations (which is simply < M | Q1 Q2 | M >
        # which Q1 and Q2 being simple deltas at LO.
        color_connections = {'LO': [
                {'open_indices': [],
                 # Information related to the first color current string applied
                 'color_string_Q1'      : color_algebra.ColorString(),
                 # The replace indices dictionary acts like the delta in color space.
                 # It is more efficient to do it like this than with actual color factors 
                 'replace_indices_Q1'   : {},
                 # Keeps track of the last dummy index used when building Q1
                 'last_dummy_index_Q1'  : self.Q1_DUMMY_STARTING_INDEX,
                 # Information related to the second color current string applied
                 'color_string_Q2'      : color_algebra.ColorString(),
                 # The replace indices dictionary acts like the delta in color space.
                 # It is more efficient to do it like this than with actual color factors
                 'replace_indices_Q2'   : {},
                 # Keeps track of the last dummy index used when building Q1
                 'last_dummy_index_Q2'  : self.Q2_DUMMY_STARTING_INDEX,
                 # The tuple representation (common to Q1 and Q2) as specified in the doc
                 # of the function 'generate_all_color_connections'
                 'tuple_representation' : (),
                 # This integer corresponds to the number of different ways in which one can
                 # arrive at this connection. For example, the tuple_representation:
                 #    ((4,-2,4),(3,-1,3))
                 # could also have been written
                 #    ((3,-2,3),(4,-1,4))
                 # but only the former will be retained as a result of the ordering principle
                 # of ascending first indices in all the triplets part of the tuple representation
                 # Therefore, n_representative will be computed to be 2 for the above connection
                 'n_representatives'    : 1}
            ]
        }

        # Each open index is stored as a 2-tuple with the following format:
        #    ( <open_index_number>, <SU(3) representation, taking values in [3,-3,8]> )
        # The operators entry corresponds to the list of color object instances building the correlator
        # and the tuple representation corresponds to the convention described in the doc of this function.
        
        for leg in process_legs:
            # Obtain the color representation of the index carried by this leg
            color_charge = model.get_particle(leg.get('id')).get_color()
            leg_number = leg.get('number')
            if color_charge == 1:
                continue
            # Initial state anti-quark and final state quarks both carry anti-fundamental color indices
            # The anti-charge generator takes a minus sign.
            elif color_charge == -3 and not leg.get('state') or color_charge == 3 and leg.get('state'):
                color_connections['LO'][0]['open_indices'].append((leg_number,-3))
            # Initial state quark and final state anti-quarks both carry fundamental color indices
            elif color_charge == 3 and not leg.get('state') or color_charge == -3 and leg.get('state'):
                color_connections['LO'][0]['open_indices'].append((leg_number,3))
            # For gluon self-interactions, we chose the second index 'b' of f^{abc} to be the one carrying 
            # "emitted" gluon's color.
            # The dual representation of the color charge takes an imaginary factor 'i', this is also what insures that
            # the color matrix remains real for the correlator of the type f * T
            elif color_charge == 8:
                color_connections['LO'][0]['open_indices'].append((leg_number,8))
            # In all cases for color-charged particle we start with deltas
            color_connections['LO'][0]['replace_indices_Q1'][leg_number]\
                = self.SINGLE_PRIME_POSITIVE_INDICES_OFFSET+leg_number
            color_connections['LO'][0]['replace_indices_Q2'][self.SINGLE_PRIME_POSITIVE_INDICES_OFFSET+leg_number]\
                = self.DOUBLE_PRIME_POSITIVE_INDICES_OFFSET+leg_number

        def create_connection_copy(correlator):
            copy =   { 'open_indices'           : list(correlator['open_indices']),
                       'color_string_Q1'        : correlator['color_string_Q1'].create_copy(),
                       'replace_indices_Q1'     : dict(correlator['replace_indices_Q1']),
                       'last_dummy_index_Q1'    : correlator['last_dummy_index_Q1'],
                       'color_string_Q2'        : correlator['color_string_Q2'].create_copy(),
                       'replace_indices_Q2'     : dict(correlator['replace_indices_Q2']),
                       'last_dummy_index_Q2'    : correlator['last_dummy_index_Q2'],
                       'tuple_representation'   : list(correlator['tuple_representation']),
                       'n_representatives'      : correlator['n_representatives']
                     }
            if 'emitted_numbers_pool' in correlator:
                copy['emitted_numbers_pool'] = list(correlator['emitted_numbers_pool'])
            return copy
        
        # We can now add the correlators for each perturbative order up to max_order
        for curr_order in range(1,max_order+1):
            
            # Now iteratively generate all n emissions, with n=curr_order. Everytime storing
            # the intermediate list of correlators generated after the i^th emission in
            # intermediate_step, which is initialized to the a lsit of a single element:
            # the LO correlator so as to have access to the list of open indices
            connections_intermediate_list = [create_connection_copy(color_connections['LO'][0])]
            # We will add an additional helper entry in the intermediate_step records which
            # is the list of emitted numbers that must still be used (initialized to [-1,-2]
            # for NNLO for instance)
            connections_intermediate_list[0]['emitted_numbers_pool'] = range(-1,-curr_order-1,-1)
            for i_emission in range(curr_order):
                next_connections_list = []
                # Append an emission to all previously generated connections
                for connection in connections_intermediate_list:
                    # Choose one open index from which to emit
                    for emitting_index, emitting_repr in connection['open_indices']:
                        # Choose one emitted index in the pool of available ones
                        for emitted_index in connection['emitted_numbers_pool']:
                            
                            # Now append the corresponding emission operator and tuple 
                            # representation to the connection from which we emit
                            new_connection = create_connection_copy(connection)
                            
                            self.add_splitting_to_connection(new_connection,
                                 emitting_index, emitting_repr, emitted_index, qqbar=False)
                            next_connections_list.append(new_connection)
                            
                            # If the emitting index is unresolved (i.e. negative), 
                            # we must allow the g > q q~ splitting. 
                            if emitting_repr == 8 and emitting_index < 0:
                                qqbar_connection = create_connection_copy(connection)
                                self.add_splitting_to_connection(qqbar_connection,
                                  emitting_index, emitting_repr, emitted_index, qqbar=True)
                                next_connections_list.append(qqbar_connection)             
                                
                # Now update the intermediate connections_intermediate list to be used for 
                # the next emission. Filter here the connections generated by ordering
                # the tuple representation and keeping only unique occurences.
                filtered_connections = {}
                for connection in next_connections_list:
                    
                    connection['tuple_representation'] = tuple(sorted(
                        connection['tuple_representation'], key = lambda _: _[0], reverse=True
                        ))
                    # Comment the line above and uncomment the next three lines below if you don't want
                    # to simplify the list of color connections using the ascending ordering
                    # principle of the first index
                    ##connection['tuple_representation'] = tuple(connection['tuple_representation'])
                    ##filtered_connections[connection['tuple_representation']] = connection
                    ##continue
                    
                    if connection['tuple_representation'] not in filtered_connections:
                        # If it wasn't already listed, simply attach it to the list of connections from
                        # which to perform the next emission
                        filtered_connections[connection['tuple_representation']] = connection
                    else:
                        # If it was already listed, simply add the number of copies this connection
                        # held to the number of copies associated to the matching connection.
                        filtered_connections[connection['tuple_representation']]['n_representatives'] += connection['n_representatives']
                
                connections_intermediate_list = filtered_connections.values()

            # Remove the temporary entry 'emitted_numbers_pool' to the generated connections
            # and register them in the main dictionary color_connections
            for corr in connections_intermediate_list:
                assert(len(corr['emitted_numbers_pool'])==0)
                del corr['emitted_numbers_pool']
            # Add this list of correlators, sorted according to their tuple representation.
            color_connections['N'*curr_order+'LO'] = sorted(connections_intermediate_list,
                key = lambda el: el['tuple_representation'])
        
        return color_connections

    def set_entries_for_correlator(self, cache, all_colored_indices_replacement,
                 color_connection_Q1, color_connection_Q2, Nc, Nc_power_min, Nc_power_max):
        """ Given the correlator, compute the entry of this color matrix."""
        
        # Useful shorthand to instantiate a ColorString from an immutable representation.
        def from_immutable(CB):
            a_cs = color_algebra.ColorString()
            a_cs.from_immutable(CB)
            return a_cs
        
        debug = False
        
        ##if ( color_connection_Q1['tuple_representation'] == ((4,-2,4 ),(4,-1,4)) and \
        ##     color_connection_Q2['tuple_representation'] == ((4,-1,4 ),(3,-2,3)) ):
        ##    print('Doing 12 vs 15: %s vs %s'%(color_connection_Q1['tuple_representation'],color_connection_Q2['tuple_representation']))
        ##    debug = True
        ##if ( color_connection_Q1['tuple_representation'] == ((4,-2,4 ),(3,-1,3)) and \
        ##     color_connection_Q2['tuple_representation'] == ((4,-1,4 ),(4,-2,4)) ):
        ##    print('Doing 11 vs 16: %s vs %s'%(color_connection_Q1['tuple_representation'],color_connection_Q2['tuple_representation']))
        ##    debug = True
        
        if debug: print("\n>>>>>Now doing '%s' vs '%s'"%(
            str(color_connection_Q1['tuple_representation']),
            str(color_connection_Q2['tuple_representation']),
            ))
        if debug: print('color_string_Q1=',color_connection_Q1['color_string_Q1'])
        if debug: print('color_string_Q2=',color_connection_Q2['color_string_Q2'])
        if debug: print('replace_indices_Q1=',color_connection_Q1['replace_indices_Q1'])
        if debug: print('replace_indices_Q2=',color_connection_Q2['replace_indices_Q2'])

        # Check if the two correlators can be squared, it might not be the case
        # when in presence of g > q q~ splittings
        if debug:print('open_indices_Q1=',sorted(color_connection_Q1['open_indices']))
        if debug:print('open_indices_Q2=',sorted(color_connection_Q2['open_indices']))
        # The radiated q and q~ from Q2 must be considered with opposite quantum numbers
        # because they are linked to those of Q1 which are not complex conjugated.
        open_indices_Q2 = sorted((ind, -1*repr if ind<0 and abs(repr)==3 else repr) 
                                     for ind, repr in color_connection_Q2['open_indices'])
        if sorted(color_connection_Q1['open_indices']) != open_indices_Q2:
            if debug: print('Incompatible quantum numbers for this correlator.')
            # The two representations are incompatible, set the color matrix to zero
            self.col_matrix_fixed_Nc = None
            self.inverted_col_matrix = None
            for i in range(len(self._col_basis1)):
                for j in range(len(self._col_basis2)):
                    self[(i, j)] = None
            return

        # Loop over all entries and compute them
        for i, CB_left in enumerate(sorted(self._col_basis1.keys())):
            for j, CB_right in enumerate(sorted(self._col_basis2.keys())):
                if self.is_symmetric and j < i:
                    continue
                if debug: print('='*80)
                # First fix negative indices that could be repeated on both sides
                CB_right = ColorMatrix.fix_summed_indices(CB_left, CB_right)
                # Convert the immutable representations to color strings, and apply
                # complex conjugation.
                if debug: print('CB_right=',CB_right)
                CB_right_CS = from_immutable(CB_right).complex_conjugate()
                if debug: print('CB_right_CS=',CB_right_CS)
                # Replace the color indices of the conjugated matrix elements by their
                # double prime image
                CB_right_CS.replace_indices( all_colored_indices_replacement )
                # We can now multiply the left Color Structure by the first connection and
                # replace all indices in CB_left_CS that are not part of the color_connection_Q1
                CB_left_CS = from_immutable(CB_left)
                if debug: print('CB_left_CS=',CB_left_CS)
                CB_left_CS.replace_indices(color_connection_Q1['replace_indices_Q1'])
                # We must create a local copy of color_string_Q1, otherwise the next 
                # replacement of indices will induce border effects.
                CB_left_CS.product(color_connection_Q1['color_string_Q1'].create_copy())
                # And similarly for Q2 now
                CB_left_CS.replace_indices(color_connection_Q2['replace_indices_Q2'])
                CB_left_CS.product(color_connection_Q2['color_string_Q2'].create_copy())

                # Convert the immutable representations to color strings
                final_color_string = color_algebra.ColorString()
                final_color_string.product(CB_left_CS)
                final_color_string.product(CB_right_CS)
                
                # Build a canonical representation of the computation to be carried in the
                # hope of recycling its result later
                canonical_entry, _ = final_color_string.to_canonical()
                if debug: print('canonical_entry=',canonical_entry,' with repl=',_)
                try:
                    # If this has already been calculated, use the result
                    result, result_fixed_Nc = cache[canonical_entry]
                    # misc.sprint('<%d| %s | %s |%d> = %s = %s'%(
                    #   i, color_connection_Q1['tuple_representation'], 
                    #   color_connection_Q2['tuple_representation'], j, 'recycled!', res_fixed_Nc))
                except KeyError:                 
                    
                    # Now simplify and evaluate the corresponding color chain
                    col_fact = color_algebra.ColorFactor([final_color_string])
                    if debug: print(col_fact)                        
                    result = col_fact.full_simplify()
                    if debug: print('result:',result)

                    # Keep only terms with Nc_max >= Nc power >= Nc_min
                    if Nc_power_min is not None:
                        result[:] = [col_str for col_str in result if col_str.Nc_power >= Nc_power_min]
                    if Nc_power_max is not None:
                        result[:] = [col_str for col_str in result if col_str.Nc_power <= Nc_power_max]
                    
                    # Set Nc to a numerical value
                    result_fixed_Nc = result.set_Nc(Nc)
                    if debug: print('result_fixed_Nc:',result_fixed_Nc)

                    if result_fixed_Nc[1] != 0:
                        raise MadGraph5Error("The elements of the color correlated matrices should always be real."+
                          " It turned out not to be the case when considering the correlator [ %s | %s ], giving %s"%
                            (color_connection_Q1['tuple_representation'], color_connection_Q2['tuple_representation'],
                                                                     str(result_fixed_Nc)))
                   
                    # Store result
                    cache[canonical_entry] = (result, result_fixed_Nc)
                    # misc.sprint('<%d| %s | %s |%d> = %s = %s'%(
                    #   i, color_connection_Q1['tuple_representation'], 
                    #   color_connection_Q2['tuple_representation'], j, 'computed!', res_fixed_Nc))
                
                # Now that we have recovered our result_fixed_Nc, we can store it in the matrix
                # Store the full result.
                self[(i, j)] = result
                if self.is_symmetric:
                    self[(j, i)] = result
                    
                # the fixed Nc one.
                self.col_matrix_fixed_Nc[(i, j)] = result_fixed_Nc
                if self.is_symmetric:
                    self.col_matrix_fixed_Nc[(j, i)] = result_fixed_Nc

                # and update the inverted dict
                if result_fixed_Nc in self.inverted_col_matrix:
                    self.inverted_col_matrix[result_fixed_Nc].append((i,j))
                    if self.is_symmetric:
                        self.inverted_col_matrix[result_fixed_Nc].append((j,i))
                else:
                    self.inverted_col_matrix[result_fixed_Nc] = [(i, j)]
                    if self.is_symmetric:
                        self.inverted_col_matrix[result_fixed_Nc] = [(j, i)]


    def build_color_correlated_matrices(self, process_legs, model, order='NLO', 
                                                Nc=3, Nc_power_min=None, Nc_power_max=None):
        """ Computes the color matrices for all relevant color connections at a given order.
        This function needs to know the process legs and a model instance so as to retrieve
        their color charges.
        For 'NLO', the returned value is a dictionary with, as keys, the tuple
           (connected_leg_number_one, connected_leg_number_one)
        and as values the tuple:
           (canonical_color_connection_representation, color_matrix)
        where the color matrix is an instance of ColorMatrix.
        """

        assert (len(process_legs)<10), "Sillyness, will mess up with color indices threshold."

        # Now list and identify all relevant color connections, which is a dictionary
        # whose format is explained in the documentation of generate_all_color_connections
        color_connections = self.generate_all_color_connections(process_legs, model, order=order)

        # All color-correlated color matrices
        all_color_correlated_matrices = {}
        # Cache results
        canonical_dict = {}
        
        # Identify all colored indices to be replaced by their double-prime image in the
        # conjugated matrix element
        all_colored_indices_replacement = {}
        for leg in process_legs:
            if model.get_particle(leg.get('id')).get_color() != 1:
                all_colored_indices_replacement[leg.get('number')] = \
                                self.DOUBLE_PRIME_POSITIVE_INDICES_OFFSET+leg.get('number')

        connection_index_offset = 1
        for curr_order in range(1,order.count('N')+1):
            order_key = 'N'*curr_order+'LO'
            
            # Loop over all possible tuples (connection_Q1, connection_Q2)
            for Q1_index, connection_Q1 in enumerate(color_connections[order_key]):
                for Q2_index, connection_Q2 in enumerate(color_connections[order_key]):
                    # Unlike at LO, there is no symmetry between (Q1, Q2) and (Q2, Q1) anymore!
                    #if Q2_index < Q1_index:
                    #    continue
                    color_correlator_identifier = (
                        connection_index_offset + Q1_index,
                        connection_index_offset + Q2_index )
                    
                    # Instantiate the ColorMatrix that will correspond to that color correlator
                    color_matrix = ColorMatrix(self._col_basis1,  
                        col_basis2 = self._col_basis2 if not self.is_symmetric else None, 
                        Nc = Nc, Nc_power_min = Nc_power_min, Nc_power_max = Nc_power_max,
                        automatic_build = False)

                    # This follows very closely what "build_matrix" does but it is unfortunately sufficiently
                    # different that it is to cumbersome to modify the color matrix above so as to reuse the build_matrix() function.
                    # The idea would be to modify the basis elements of these color matrix so as to include there directly
                    # the color connection, but it has many pitfalls, one of which being that it could alter the 'sorted' order
                    # of the color basis keys, on which we rely for the definition of the color matrix ported in the ME code.
                    color_matrix.set_entries_for_correlator( canonical_dict, 
                        all_colored_indices_replacement, connection_Q1, connection_Q2, 
                                               Nc=3, Nc_power_min=None, Nc_power_max=None )

                    # And we can now add our finalized color_matrix to the dictionary that will be returned
                    all_color_correlated_matrices[color_correlator_identifier] = ( 
                        ( str(connection_Q1['color_string_Q1']), 
                          str(connection_Q2['color_string_Q2']) ), color_matrix )
            
            # Increase the offset by the length of the list of connections up to this point
            connection_index_offset += len(color_connections[order_key])
            
        return all_color_correlated_matrices, color_connections

    def build_matrix(self, Nc=3,
                     Nc_power_min=None,
                     Nc_power_max=None,
                     is_symmetric=False):
        """Create the matrix using internal color basis objects. Use the stored
        color basis objects and takes Nc and Nc_min/max parameters as __init__.
        If is_isymmetric is True, build only half of the matrix which is assumed
        to be symmetric."""

        canonical_dict = {}
        
        for i1, struct1 in \
                    enumerate(sorted(self._col_basis1.keys())):
            for i2, struct2 in \
                    enumerate(sorted(self._col_basis2.keys())):
                # Only scan upper right triangle if symmetric
                if is_symmetric and i2 < i1:
                    continue

                # Fix indices in struct2 knowing summed indices in struct1
                # to avoid duplicates
                new_struct2 = self.fix_summed_indices(struct1, struct2)

                # Build a canonical representation of the two immutable struct
                canonical_entry, dummy = \
                            color_algebra.ColorString().to_canonical(struct1 + \
                                                                   new_struct2)

                try:
                    # If this has already been calculated, use the result
                    result, result_fixed_Nc = canonical_dict[canonical_entry]
                except KeyError:
                    # Otherwise calculate the result
                    result, result_fixed_Nc = \
                            self.create_new_entry(struct1,
                                                  new_struct2,
                                                  Nc_power_min,
                                                  Nc_power_max,
                                                  Nc)
                    # Store both results
                    canonical_dict[canonical_entry] = (result, result_fixed_Nc)

                # Store the full result...
                self[(i1, i2)] = result
                if is_symmetric:
                    self[(i2, i1)] = result

                # the fixed Nc one ...
                self.col_matrix_fixed_Nc[(i1, i2)] = result_fixed_Nc
                if is_symmetric:
                    self.col_matrix_fixed_Nc[(i2, i1)] = result_fixed_Nc
                # and update the inverted dict
                if result_fixed_Nc in self.inverted_col_matrix.keys():
                    self.inverted_col_matrix[result_fixed_Nc].append((i1,
                                                                      i2))
                    if is_symmetric:
                        self.inverted_col_matrix[result_fixed_Nc].append((i2,
                                                                          i1))
                else:
                    self.inverted_col_matrix[result_fixed_Nc] = [(i1, i2)]
                    if is_symmetric:
                        self.inverted_col_matrix[result_fixed_Nc] = [(i2, i1)]

    def create_new_entry(self, struct1, struct2,
                         Nc_power_min, Nc_power_max, Nc):
        """ Create a new product result, and result with fixed Nc for two color
        basis entries. Implement Nc power limits."""

        # Create color string objects corresponding to color basis 
        # keys
        col_str = color_algebra.ColorString()
        col_str.from_immutable(struct1)

        col_str2 = color_algebra.ColorString()
        col_str2.from_immutable(struct2)

        # Complex conjugate the second one and multiply the two
        col_str.product(col_str2.complex_conjugate())

        # Create a color factor to store the result and simplify it
        # taking into account the limit on Nc
        col_fact = color_algebra.ColorFactor([col_str])
        result = col_fact.full_simplify()

        # Keep only terms with Nc_max >= Nc power >= Nc_min
        if Nc_power_min is not None:
            result[:] = [col_str for col_str in result \
                         if col_str.Nc_power >= Nc_power_min]
        if Nc_power_max is not None:
            result[:] = [col_str for col_str in result \
                         if col_str.Nc_power <= Nc_power_max]

        # Calculate the fixed Nc representation
        result_fixed_Nc = result.set_Nc(Nc)

        return result, result_fixed_Nc

    def __str__(self):
        """Returns a nicely formatted string with the fixed Nc representation
        of the current matrix (only the real part)"""

        mystr = '\n\t' + '\t'.join([str(i) for i in \
                                    range(len(self._col_basis2))])

        for i1 in range(len(self._col_basis1)):
            mystr = mystr + '\n' + str(i1) + '\t'
            mystr = mystr + '\t'.join(['%i/%i' % \
                        (self.col_matrix_fixed_Nc[(i1, i2)][0].numerator,
                        self.col_matrix_fixed_Nc[(i1, i2)][0].denominator) \
                        for i2 in range(len(self._col_basis2))])

        return mystr

    def get_line_denominators(self):
        """Get a list with the denominators for the different lines in
        the color matrix"""
        
        # It can happen that the color matrix for fixed NC is set to None, in cases
        # of color connections with vanishing interferences when computing color-correlated
        # matrix elements (for higher order computations). The convention in this case
        # is to return denominators which are all zeros
        if self.col_matrix_fixed_Nc is None:
            return [0 for i1 in range(len(self._col_basis1))]
        else:
            den_list = []
            for i1 in range(len(self._col_basis1)):
                den_list.append(self.lcmm(*[\
                            self.col_matrix_fixed_Nc[(i1, i2)][0].denominator for \
                                            i2 in range(len(self._col_basis2))]))
        return den_list

    def get_line_numerators(self, line_index, den):
        """Returns a list of numerator for line line_index, assuming a common
        denominator den."""

        # It can happen that the color matrix for fixed NC is set to None, in cases
        # of color connections with vanishing interferences when computing color-correlated
        # matrix elements (for higher order computations). The convention in this case
        # is to return all zero numerators
        if self.col_matrix_fixed_Nc is None:
            return [0 for i2 in range(len(self._col_basis2))]
        else:
            return [self.col_matrix_fixed_Nc[(line_index, i2)][0].numerator * \
                den / self.col_matrix_fixed_Nc[(line_index, i2)][0].denominator \
                for i2 in range(len(self._col_basis2))]

    @classmethod
    def fix_summed_indices(self, struct1, struct2):
        """Returns a copy of the immutable Color String representation struct2 
        where summed indices are modified to avoid duplicates with those
        appearing in struct1. Assumes internal summed indices are negative."""

        # First, determines what is the smallest index appearing in struct1
        #list2 = reduce(operator.add,[list(elem[1]) for elem in struct1])
        list2 = sum((list(elem[1]) for elem in struct1),[])
        if not list2: 
            min_index = -1
        else:
           min_index = min(list2) - 1

        # Second, determines the summed indices in struct2 and create a 
        # replacement dictionary
        repl_dict = {}
        #list2 = reduce(operator.add,
        #               [list(elem[1]) for elem in struct1])
        for summed_index in list(set([i for i in list2 \
                                      if list2.count(i) == 2])):
            repl_dict[summed_index] = min_index
            min_index -= 1

        # Three, create a new immutable struct by doing replacements in struct2
        return_list = []
        for elem in struct2:
            fix_elem = [elem[0], []]
            for index in elem[1]:
                try:
                    fix_elem[1].append(repl_dict[index])
                except Exception:
                    fix_elem[1].append(index)
            return_list.append((elem[0], tuple(fix_elem[1])))

        return tuple(return_list)

    @staticmethod
    def lcm(a, b):
        """Return lowest common multiple."""
        return a * b // fractions.gcd(a, b)

    @staticmethod
    def lcmm(*args):
        """Return lcm of args."""
        if args:
            return reduce(ColorMatrix.lcm, args)
        else:
            return 1

