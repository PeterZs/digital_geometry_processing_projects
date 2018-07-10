#!/bin/sh

EXEC=./bin/remesh.exe.opt
make remesh METHOD=opt


# # --------------------- Sphere
# $EXEC -i meshes/report/sphere.off -o results/data/sphere_b r_obtuse=1 
# $EXEC -i meshes/report/sphere.off -o results/data/sphere_c n_verts=10000

# $EXEC -i meshes/report/sphere.off -o results/data/sphere_d proj_type=current
# $EXEC -i meshes/report/sphere.off -o results/data/sphere_e 

# # --------------------- Cow
# $EXEC -i meshes/report/cow1.obj -o results/data/cow_b n_verts=6000
# $EXEC -i meshes/report/cow1.obj -o results/data/cow_c n_verts=20000

# # --------------------- Hand
# $EXEC -i meshes/report/hand.off -o results/data/hand_b r_obtuse=1
# $EXEC -i meshes/report/hand.off -o results/data/hand_c n_verts=20000

# # --------------------- Face
# $EXEC -i meshes/report/nicolo.off -o results/data/nicolo_b r_obtuse=1
# $EXEC -i meshes/report/nicolo.off -o results/data/nicolo_c n_verts=20000

# --------------------- Horse
# $EXEC -i meshes/report/horse.obj -o results/data/horse_b n_verts=10000
# $EXEC -i meshes/report/horse.obj -o results/data/horse_c n_verts=5000
# $EXEC -i meshes/report/horse.obj -o results/data/horse_e n_verts=10000 proj_type=current
# $EXEC -i meshes/report/horse.obj -o results/data/horse_f n_verts=5000 proj_type=current


./bin/remesh.exe.opt n_verts=5000 -i meshes/maxplanck_s_m1_30000.off -o results/data/mp_a #4568
./bin/remesh.exe.opt n_verts=10000 -i meshes/maxplanck_s_m1_30000.off -o results/data/mp_b #4210
./bin/remesh.exe.opt n_verts=15000 -i meshes/maxplanck_s_m1_30000.off -o results/data/mp_c #3894
./bin/remesh.exe.opt n_verts=20000 -i meshes/maxplanck_s_m1_30000.off -o results/data/mp_d # 3542
./bin/remesh.exe.opt n_verts=30000 -i meshes/maxplanck_s_m1_30000.off -o results/data/mp_e #3520
