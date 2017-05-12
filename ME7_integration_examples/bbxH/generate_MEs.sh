rm -rf bbxH_tree_MEs
rm -rf bbxH_loop_ME
cp restrict_no_b_mass_with_yukawa.dat ../../models/loop_sm/restrict_no_b_mass_with_yukawa.dat
../../bin/mg5_aMC ME_generation.mg5
mv bbxH_loop_ME/Cards bbxH_loop_ME/orig_Cards
cd bbxH_loop_ME; ln -s ../Cards Cards; cd -;
./compile_MEs.sh
rm -f nsqso_born.inc
rm -f additional_command
