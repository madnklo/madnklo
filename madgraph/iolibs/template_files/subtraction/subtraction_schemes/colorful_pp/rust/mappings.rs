    @classmethod
    def map_to_lower_multiplicity(
        cls, PS_point, singular_structure, momenta_dict, squared_masses=None, kinematic_variables=None, compute_jacobian=False ):
        """Map the resolved momenta of a real emission process to on-shell kinematics by recoling against the initial state partons.

        :param PS_point: dictionary of Real kinematics
        :param singular_structure: SingularStructure describing the limit we want to factorize with this mapping
        :param momenta_dict: Routing of momenta label
        :param squared_masses: dictionary
        :param kinematic_variables: dictionary to store the variables of the unresolved PS
        :param compute_jacobian: Boolean, specify that the jacobian should be computed if it is necessary to divide by it in the local CT.
        :return:
            * new_PS_point: a dictionary with Born kinematics
            * mapping_variables: a dictionary with variables, possible used in the current evaluation.
        """

        # Consistency checks
        assert isinstance(momenta_dict, sub.bidict)
        if not cls.is_valid_structure(singular_structure):
            raise MappingError("Singular structure '%s' is not supported by mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Recoilers must be the two initial state legs:
        if len(singular_structure.legs)!=2 and not all(leg.state==leg.INITIAL for leg in singular_structure.legs):
            raise MappingError("Singular structure '%s' legs should be two initial states when using the mapping '%s'"%(
                                                    str(singular_structure), cls.__name__))

        # Construct a new phase space point
        new_PS_point = PS_point.get_copy()

        # Build the total collinear momentum and remove all collinear legs from the mapped PS point
        # Below must be upgraded for nested collinears
        pC = LorentzVector()
        collinear_children = set([])
        for substructure in singular_structure.substructures: # S(3,4).legs returns [3,4]
            children = tuple(leg.n for leg in substructure.legs)
            for child in children:
                pC += new_PS_point.pop(child)
            collinear_children |= set(children)

        # Build the total momentum Q from the initial states. This is of course only applicable for 2>N kinematics.
        # Fetch the initial states numbers using the recoiler ids.
        Q = LorentzVector()
        Q += PS_point[singular_structure.legs[0].n]
        Q += PS_point[singular_structure.legs[1].n]
        Q_square = Q.square()
        pa_dot_pb = Q_square/2.

        # Compute the alpha parameter
        alpha = 0.5*( pC.dot(Q)/pa_dot_pb - math.sqrt( (pC.dot(Q)/pa_dot_pb)**2 - 4.* (pC.square()/Q_square) ) )

        # Recoil against initial states
        for recoiler in singular_structure.legs:
            new_PS_point[recoiler.n] = PS_point[recoiler.n]*(1.-alpha)

        # Find the leg mother label
        # Note: in the case of nested collinear limits, the fetching below must be upgraded
        mother_label = momenta_dict.inv[frozenset(collinear_children)]
        new_PS_point[mother_label] = pC - alpha*Q

        mapping_variables = {'Q': Q, 'alpha': alpha}

        # The jacobian may be one, but not sure; haven't computed it yet.
        if compute_jacobian:
            raise NotImplementedError

        return new_PS_point, mapping_variables
