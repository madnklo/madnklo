struct Variables {
    zs: Vec<f64>,
    kTs: Vec<LorentzVector<f64>>,
    ss: Vec<Vec<f64>>,
    ss_i_others: Vec<Vec<f64>>,
}

/// Variables for 'n' *pure final state* collinear recoiling exclusively against the initial state.
fn colorful_pp_FFn_variables(higher_PS_point: HashMap<usize, LorentzVector<f64>>, 
    qC: LorentzVector<f64>, children: &[usize], Q: LorentzVector<f64>, alpha: f64) -> Variables {

    let all_p_fs: Vec<_> = children.iter().map(|child| &higher_PS_point[*child]).collect();

    let pC = all_p_fs.iter().sum();
    let Q_square = Q.square();
    let pa_dot_pb = Q_square/2.;
    let alpha = 0.5*( pC.dot(Q)/pa_dot_pb - (pC.dot(Q)/pa_dot_pb)**2 - 4.* (pC.square()/Q_square) ).sqrt() ) ;

    // Loop over all final state momenta to obtain the corresponding variables
    let mut zs = Vec::with_capacity(all_p_fs.len());
    for p_fs in all_p_fs {
        zs.push(p_fs.dot(Q));
    }

    let normalisation = zs.sum();
    for zsi in &mut zs {
        *zsi /= normalisation;
    }

    // Add additional ss's
    let mut ss_i_j = vec![vec![0.; all_p_fs.len()]; all_p_fs.len()];
    let mut ss_i_j_tot = 0.;
    for i_fs, p_fs_i in all_p_fs.iter().enumerate() {
        for j_fs in i_fs+1..all_p_fs.len() {
            ss_i_j[i_fs][j_fs] = 2.*p_fs_i.dot(all_p_fs[j_fs]);
            ss_i_j_tot += ss_i_j[i_fs][j_fs];
            ss_i_j[j_fs][i_fs] = ss_i_j[i_fs][j_fs];
        }
    }

    let ss_i_others = all_p_fs.iter().map(|p| 2*(pC - p).dot(p)).collect();

    let xis = zs.iter().zip(ss_i_others).map(|z, so| (z-so/(alpha*(2.*pC.dot(&Q))))).collect();

    if len(all_p_fs)==3 {
        // WARNING function below is *not* symmetric in r and s.
        def bigZ(i,r,s):
            return (ss_i_j[(i,r)] - ss_i_j[(r,s)] - 2.*zs[r]*ss_i_j_tot)/(alpha*2.*qC.dot(Q))
    }
    else len(all_p_fs)==2 {
        def bigZ(i,r):
            return (ss_i_j[(i,r)]*(zs[r]-zs[i]))/(alpha*2.*qC.dot(&Q))
    }
    else:
        notimplemented()

    kTs = Vec::with_capacity(all_p_fs.len());
    for (i, p_fs_i) in all_p_fs.iter().enumerate() {
        let mut kT = LorentzVector::default();

        for j, p_fs_j in all_p_fs.iter().enumerate() {
            kT += xis[i]*p_fs_j - xis[j] * p_fs_i;
        }

        match all_p_fs.len() {
            2 => {
                let r = 1 if i == 0 else 0;
                kT += (ss_i_j[i][r]*(zs[r]-zs[i]))/(alpha*2.*qC.dot(&Q)) * qC;
            }
            3 => {
                // WARNING function below is *not* symmetric in r and s.
                let r = 1 if i == 0 else 0;
                let s = 1 if i == 2 else 2;
                let Z = (ss_i_j[i][r] - ss_i_j[r][s] - 2.*zs[r]*ss_i_j_tot)/(alpha*2.*qC.dot(&Q));
                kT += Z * qC;
            }
            x => notimplemented!("Case {} is not supported", x)
        }

        kTs[i] = kT;
    }

    Variables {
        zs,
        kTs,
        ss,
        ss_i_others
    }
}