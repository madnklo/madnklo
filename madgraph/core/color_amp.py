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

    def get_color_generator_for_leg(self, leg, model, this_leg_index_offset, other_leg_index_offset):
        """ Builds the color generator appropriate for a given leg, given its color charge.
        index offsets specify what index offset to use to differentiate the indices that
        are part of the color connections and also part of <M| and |M>.
        The leg_index_offsets argument specifies which offset to apply to select exactly 
        how they will be contracted (i.e. with the indices of the <M|, or |M> or even of the color 
        connector itself. """

        # Make sure indices summation in this color connection doesn't conflict with 
        # the ones of other ColorStrings that will be part of the same chain.
        negative_offset = -10000

        color_charge = model.get_particle(leg.get('id')).get_color()
        leg_number = leg.get('number')
        this_leg_index = this_leg_index_offset + leg_number
        other_leg_index = other_leg_index_offset + leg_number
        
        if color_charge==1:
            return None
        if color_charge not in [3,-3,8]:
            raise MadGraph5Error("Color connected matrix elements between particles color-charged"+
                                 " in the %d representation are not implemented."%abs(color_charge))

        # Initial state anti-quark and final state quarks both carry anti-fundamental color indices
        # The anti-charge generator takes a minus sign.
        if color_charge == -3 and not leg.get('state') or color_charge == 3 and leg.get('state'):
            return color_algebra.ColorString([color_algebra.T(negative_offset, other_leg_index, this_leg_index)],
                                                                                    coeff=fractions.Fraction(-1, 1))
        # Initial state quark and final state anti-quarks both carry fundamental color indices
        if color_charge == 3 and not leg.get('state') or color_charge == -3 and leg.get('state'):
            return color_algebra.ColorString([color_algebra.T(negative_offset, this_leg_index, other_leg_index)])
        # For gluon self-interactions, we chose the second index 'b' of f^{abc} to be the one carrying 
        # "emitted" gluon's color.
        # The dual representation of the color charge takes an imaginary factor 'i', this is also what insures that
        # the color matrix remains real for the correlator of the type f * T
        if color_charge == 8:
#ORIG            return color_algebra.ColorString([color_algebra.f(other_leg_index, negative_offset,this_leg_index)],
#ORIG                                                                                                 is_imaginary=True)
            return color_algebra.ColorString([color_algebra.f(this_leg_index, negative_offset, other_leg_index)],
                                                                                                 is_imaginary=True)

    def generate_all_color_connections(self, index_offset, process_legs, model, order='NLO'):
        """Returns a dictionary whose keys is the "identifier" of the color connection. At NLO, this identifier is
        a tuple showing pairs of ID's from which the particle is emitted / re-absorbed and values are the ColorString 
        corresponding to this connection.
        """

        color_connections = {}
        
        if order != 'NLO':
            raise MadGraph5Error("Color connections for expansion order '%s' not implemented yet."%order)

        for leg1 in process_legs:
            color_string_1 = self.get_color_generator_for_leg(leg1, model, 0, index_offset)
            if color_string_1 is None:
                continue
            for leg2 in process_legs:
                # Only build half of the symmetric color connection CC_ij = <M| T_i T_j |M> matrix
                # We will build the diagonal terms later
                if leg2.get('number')<=leg1.get('number'):
                    continue
#ORIG                color_string_2 = self.get_color_generator_for_leg(leg2, model, index_offset, 0)
                color_string_2 = self.get_color_generator_for_leg(leg2, model, 0, index_offset)
                if color_string_2 is None:
                    continue
                color_connection = color_string_1.create_copy()
#ORIG                color_connection.product(color_string_2.complex_conjugate())
                color_connection.product(color_string_2)
                
                color_connections[(leg1.get('number'), leg2.get('number'))] = \
                                                    (color_connection, color_connection.to_immutable())
        
        # Now build the diagonal elements which will always be simpliable in terms of delta's and the resulting color matrix
        # will be proportional to the original one, but since we want to retain the possibility of truncating Nc_orders
        # we will not take advantage of this fact here
        for leg in process_legs:
            # The indices of the color generators inserted will saturate with themselves, so we must take
            # them distinct from the ones that will appear in the matrix-element <M| and |M>, hence we use
            # 100*index_offset
            color_connection = self.get_color_generator_for_leg(leg, model, 0, 10*index_offset)
            if color_connection is None:
                continue
            
#ORIG            color_string_2 = self.get_color_generator_for_leg(leg, model, index_offset, 10*index_offset) 
#ORIG            color_connection.product(color_string_2.complex_conjugate())

            color_string_2 = self.get_color_generator_for_leg(leg, model, 10*index_offset, index_offset) 
            color_connection.product(color_string_2)            
            
#            misc.sprint(repr(color_connection), color_algebra.ColorFactor([color_connection]).full_simplify())
            color_connections[(leg.get('number'), leg.get('number'))] = \
                                                (color_connection, color_connection.to_immutable())        
        
        return color_connections

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

        # Replacement color index offset. This is because when evaluating:
        #  ColorString1(indices1) * ColorConnection(CCindices) * ColorString2(indices2)
        # The indices of indices 2 and indices1 overlap with those of CC indices and we have to 
        # differentiate them by offsetting by index_offset those in indices2 also present in 
        # CCindices and indices1.
        index_offset = 1000
        
        assert (len(process_legs)<1000), "Sillyness, will mess up with color indices threshold."

        # Useful shorthand to instantiate a ColorString from an immutable representation.
        def from_immutable(CB):
            a_cs = color_algebra.ColorString()
            a_cs.from_immutable(CB)
            return a_cs
       
        # Now list and identify all relevant color connections
        # Keys are tuple showing pairs of ID's from which the particle is emitted / re-absorbed
        # Values are the ColorString corresponding to this connection
        color_connections = self.generate_all_color_connections(
                                                index_offset, process_legs, model, order=order)
        
        # All color-correlated color matrices
        all_color_correlated_matrices = {}
        # Cache results
        canonical_dict = {}
#        misc.sprint('Total of %d color connections.'%len(color_connections))
        for color_connection_identifier in sorted(color_connections.keys()):
            # Note that this uses the "identifier" format of the color connection suited for NLO.
            # This will need to be revisited for NNLO.
            leg1_number, leg2_number = color_connection_identifier
            color_connection = color_connections[color_connection_identifier]
#            misc.sprint('Considering color connection (%d, %d) -> %s ...'%(leg1_number, leg2_number, str(color_connection[1])))
            # Instantiate the ColorMatrix that will correspond to that color connection
            color_matrix = ColorMatrix(self._col_basis1,  
                col_basis2 = self._col_basis2 if not self.is_symmetric else None, 
                Nc = Nc, Nc_power_min = Nc_power_min, Nc_power_max = Nc_power_max,
                automatic_build = False)

            # This follows very closely what "build_matrix" does but it is unfortunately sufficiently
            # different that it is to cumbersome to modifiy the color matrix above so as to reuse the build_matrix() function.
            # The idea would be to modify the basis elements of these color matrix so as to include there directly
            # the color connection, but it has many pitfalls, one of which being that it could alter the 'sorted' order
            # of the color basis keys, on which we rely for the definition of the color matrix ported in the ME code.
            for i, CB_left in enumerate(sorted(color_matrix._col_basis1.keys())):
                for j, CB_right in enumerate(sorted(color_matrix._col_basis2.keys())):
                    if color_matrix.is_symmetric and j < i:
                        continue

                    # Fix negative indices that could be repeated on both sides
                    CB_right = ColorMatrix.fix_summed_indices(CB_left, CB_right)
                    # Convert the immutable representations to color strings
                    CB_right_CS = from_immutable(CB_right).complex_conjugate()
                    CB_right_CS.replace_indices({leg1_number : index_offset + leg1_number,
                                                 leg2_number : index_offset + leg2_number})
                    # Build a canonical representation of the computation to be carried in the hope of recycling its result
                    canonical_entry, _ = color_algebra.ColorString().to_canonical(
                                                                CB_left+color_connection[1]+CB_right_CS.to_immutable())
                    try:
                        # If this has already been calculated, use the result
                        result, result_fixed_Nc = canonical_dict[canonical_entry]
                        # misc.sprint('<%d| T_%d T_%d |%d> = %s = %s'%(
                        #                  i, leg1_number, leg2_number, j, 'recycled!', res_fixed_Nc))
                    except KeyError:
                        # Convert the immutable representations to color strings
                        CB_left_CS = from_immutable(CB_left)
                        final_color_string = color_algebra.ColorString()
                        final_color_string.product(CB_left_CS)
                        final_color_string.product(color_connection[0])
                        final_color_string.product(CB_right_CS)
                        
                        
                        # Now simplify and evaluate the corresponding color chain
                        col_fact = color_algebra.ColorFactor([final_color_string])
                        result = col_fact.full_simplify()
                        # Keep only terms with Nc_max >= Nc power >= Nc_min
                        if Nc_power_min is not None:
                            result[:] = [col_str for col_str in result if col_str.Nc_power >= Nc_power_min]
                        if Nc_power_max is not None:
                            result[:] = [col_str for col_str in result if col_str.Nc_power <= Nc_power_max]
                        
                        # Set Nc to a numerical value
                        result_fixed_Nc = result.set_Nc(Nc)
                        
                        if result_fixed_Nc[1] != 0:
                            raise MadGraph5Error("The elements of the color correlated matrices should always be real."+
                                            " It turned out not to be the case when considering the correlator %s : %s"%
                                                                             (str(color_connection[1]),str(result_fixed_Nc)))
                       
                        # Store result
                        canonical_dict[canonical_entry] = (result, result_fixed_Nc)
                        # misc.sprint('<%d| T_%d T_%d |%d> = %s = %s'%(
                        #    i, leg1_number, leg2_number, j, repr(final_color_string), res_fixed_Nc))
                    
                    # Now that we have recovered our result_fixed_Nc, we can store it in the matrix
                    # Store the full result.
                    color_matrix[(i, j)] = result
                    if color_matrix.is_symmetric:
                        color_matrix[(j, i)] = result
                        
                    # the fixed Nc one.
                    color_matrix.col_matrix_fixed_Nc[(i, j)] = result_fixed_Nc
                    if color_matrix.is_symmetric:
                        color_matrix.col_matrix_fixed_Nc[(j, i)] = result_fixed_Nc

                    # and update the inverted dict
                    if result_fixed_Nc in color_matrix.inverted_col_matrix:
                        color_matrix.inverted_col_matrix[result_fixed_Nc].append((i,j))
                        if color_matrix.is_symmetric:
                            color_matrix.inverted_col_matrix[result_fixed_Nc].append((j,i))
                    else:
                        color_matrix.inverted_col_matrix[result_fixed_Nc] = [(i, j)]
                        if color_matrix.is_symmetric:
                            color_matrix.inverted_col_matrix[result_fixed_Nc] = [(j, i)]

            # And we can now add our finalized color_matrix to the dictionary that will be returned
            all_color_correlated_matrices[color_connection_identifier] = (color_connection[1], color_matrix)

        return all_color_correlated_matrices

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

        den_list = []
        for i1 in range(len(self._col_basis1)):
            den_list.append(self.lcmm(*[\
                        self.col_matrix_fixed_Nc[(i1, i2)][0].denominator for \
                                        i2 in range(len(self._col_basis2))]))
        return den_list

    def get_line_numerators(self, line_index, den):
        """Returns a list of numerator for line line_index, assuming a common
        denominator den."""

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

