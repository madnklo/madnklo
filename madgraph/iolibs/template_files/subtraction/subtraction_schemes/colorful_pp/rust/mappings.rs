use commons_current::{SingularStructure, SubtractionLegState};

///
/// Map the resolved momenta of a real emission process to on-shell kinematics by recoling against the initial state partons.
///
/// `param` PS_point: dictionary of Real kinematics
/// `param` singular_structure: SingularStructure describing the limit we want to factorize with this mapping
/// `param` momenta_dict: Routing of momenta label
/// `param` compute_jacobian: Boolean, specify that the jacobian should be computed if it is necessary to divide by it in the local CT.
/// `return`:
///     * new_PS_point: a dictionary with Born kinematics
///     * mapping_variables: a dictionary with variables, possible used in the current evaluation.
fn map_to_lower_multiplicity(
    PS_point: HashMap<usize, LorentzVector<f64>>, singular_structure: &SingularStructure, 
    momenta_dict: HashMap<Vec<usize>, usize>, compute_jacobian: bool) {

    // Recoilers must be the two initial state legs:
    if singular_structure.legs.len() !=2 && !singular_structure.legs.iter().all(|leg| leg.state== SubtractionLegState::Initial) {
        panic!("Mapping error: singular structure {:?} legs should be two initial states", singular_structure.legs);
    }

    // Construct a new phase space point
    let mut new_PS_point = PS_point.clone();

    // Build the total collinear momentum and remove all collinear legs from the mapped PS point
    // Below must be upgraded for nested collinears
    let mut pC = LorentzVector::default();
    let mut collinear_children = vec![];
    for substructure in singular_structure.substructures { // S(3,4).legs returns [3,4]
        children = tuple(leg.n for leg in substructure.legs);
        for child in substructure.legs {
            pC += new_PS_point.remove(child.leg_id).unwrap();
            if !collinear_children.contains(child.leg_id) {
                collinear_children.push(child.leg_id);
            };
        }
    }
    collinear_children.sort();

    // Build the total momentum Q from the initial states. This is of course only applicable for 2>N kinematics.
    // Fetch the initial states numbers using the recoiler ids.
    let mut Q = LorentzVector::default();
    Q += PS_point[singular_structure.legs[0].leg_id];
    Q += PS_point[singular_structure.legs[1].leg_id];
    let Q_square = Q.square();
    let pa_dot_pb = Q_square/2.;

    // Compute the alpha parameter
    let alpha = 0.5*( pC.dot(Q)/pa_dot_pb - ( (pC.dot(Q)/pa_dot_pb).powi(2) - 4.* (pC.square()/Q_square) ).sqrt() );

    // Recoil against initial states
    for recoiler in singular_structure.legs {
        new_PS_point.insert(recoiler.leg_id, PS_point[recoiler.leg_id]*(1.-alpha));
    }

    // Find the leg mother label
    // Note: in the case of nested collinear limits, the fetching below must be upgraded
    let mother_label = momenta_dict[collinear_children];
    new_PS_point.insert(mother_label, pC - alpha*Q);

    // The jacobian may be one, but not sure; haven't computed it yet.
    if compute_jacobian {
        notimplemented!("Computation of the Jacobian is not implemented");
    }

    (new_PS_point, 1.0)
}
