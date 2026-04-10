from pymol import cmd
import wobj
obj = wobj.import_wobj("output.obj")
cmd.load("../step5_input.pdb", "gromacs_step5")
cmd.remove("gromacs_step5 and !polymer")
wobj.draw_wobj(obj, groupname = "avg_density")
wobj.draw_wobj(obj, groupname = "avg_pf_hydrophobicity")
wobj.draw_wobj(obj, groupname = "avg_energy")
wobj.draw_wobj(obj, groupname = "avg_pl_hydrophobicity")
wobj.draw_wobj(obj, groupname = "avg_radius")
cmd.save('test.pse')

