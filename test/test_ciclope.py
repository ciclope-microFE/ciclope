import ciclope.core.voxelFE
from ciclope.utils.recon_utils import read_tiff_stack
from ciclope.utils.preprocess import remove_unconnected
import ciclope
import numpy as np
import unittest

test_data_file = 'trab/slice_00.tif'

class run_tests(unittest.TestCase):
    def test_remove_unconnected(self):
        # read and segment test dataset
        pippo_BW = read_tiff_stack(test_data_file) > 100
        pippo_BW_remove = remove_unconnected(pippo_BW)
        self.assertEqual((np.sum(pippo_BW) - np.sum(pippo_BW_remove)), 4)

    def test_tetraFE_shell_mesh_pymcubes(self):
        # (pymcubes method)
        vert, tri, mesh = ciclope.tetraFE.shell_mesh(remove_unconnected((read_tiff_stack(test_data_file)>100)), method='pymcubes', voxelsize=[0.12, 0.12, 0.12])
        # check if 3D shell mesh fields (vertices, triangles) were produced
        self.assertEqual(vert.shape[1], 3)
        self.assertEqual(tri.shape[1], 3)

    def test_tetraFE_shell_mesh_pygalmesh(self):
        # (pygalmesh method)
        vert, tri, mesh = ciclope.tetraFE.shell_mesh(remove_unconnected((read_tiff_stack(test_data_file)>100)), method='pygalmesh', voxelsize=[0.12, 0.12, 0.12])
        self.assertNotEqual(len(mesh.points), 0)
        self.assertEqual(mesh.cells[0][0], 'triangle')
        self.assertNotEqual(len(mesh.cells[0][1]), 0)

    def test_tetraFE_cgal_mesh(self):
        mesh = ciclope.tetraFE.cgal_mesh(remove_unconnected((read_tiff_stack(test_data_file)>100)), [0.12, 0.12, 0.12], meshtype='tetra')
        # check that points and cells fields are not empty
        self.assertNotEqual(len(mesh.points), 0)
        self.assertEqual(mesh.cells[0][0], 'tetra')
        self.assertNotEqual(len(mesh.cells[0][1]), 0)

    def test_mesh2tetrafe(self):
        ciclope.tetraFE.mesh2tetrafe(ciclope.tetraFE.cgal_mesh(remove_unconnected((read_tiff_stack(test_data_file)>100)), [0.12, 0.12, 0.12], meshtype='tetra'), 'test_input.inp', 'foo.inp')
        # read and check input CALCULIX file: check existence of essential fields
        with open("foo.inp", 'r') as f:
            lines = f.readlines()
            node = False
            element = False
            elset = False
            elset_z1 = False
            nset_z0 = False
            step = False
            bc = False

            for line in lines:
                if line.find('*NODE') != -1:
                    node = True
                if line.find('*ELEMENT, TYPE=C3D4') != -1:
                    element = True
                if line.find('*ELSET, ELSET=SETALL') != -1:
                    elset = True
                if line.find('*ELSET, ELSET=ELEMS_Z1') != -1:
                    elset_z1 = True
                if line.find('*NSET, NSET=NODES_Z0') != -1:
                    nset_z0 = True
                if line.find('*STEP') != -1:
                    step = True
                if line.find('NODES_T, 3, 3, -0.002') != -1:
                    bc = True

        self.assertEqual(node, True)
        self.assertEqual(element, True)
        self.assertEqual(elset, True)
        self.assertEqual(elset_z1, True)
        self.assertEqual(nset_z0, True)
        self.assertEqual(step, True)
        self.assertEqual(bc, True)

    def test_voxelFE_vol2ugrid(self):
        pippo_mesh = ciclope.core.voxelFE.vol2ugrid(read_tiff_stack(test_data_file), [0.12,0.12,0.12], GVmin=100)
        self.assertEqual(pippo_mesh.cells[0].type, "hexahedron")
        self.assertEqual(pippo_mesh.cells[0].data.shape, (209,8))
        self.assertEqual(pippo_mesh.cells[0].data.shape, (209,8))
        self.assertEqual(pippo_mesh.points.shape, (1331,3))
        self.assertEqual(pippo_mesh.point_sets['NODES_X0Y0Z0'], [18,19])
        self.assertEqual((pippo_mesh.points[100] == np.array([0.12, 1.08, 0])).all(), True)

    def test_mesh2voxelFE(self):
        ciclope.voxelFE.mesh2voxelfe(ciclope.core.voxelFE.vol2ugrid(read_tiff_stack(test_data_file), [0.12,0.12,0.12], GVmin=100), 'test_input.inp', 'foo.inp')
        # read and check input CALCULIX file: check existence of essential fields
        with open("foo.inp", 'r') as f:
            lines = f.readlines()
            node = False
            element = False
            elset_z1 = False
            nset_z0 = False
            step = False
            bc = False

            for line in lines:
                if line.find('*NODE') != -1:
                    node = True
                if line.find('*ELEMENT, TYPE=C3D8, ELSET=SET255') != -1:
                    element = True
                if line.find('*ELSET, ELSET=CELLS_Z1') != -1:
                    elset_z1 = True
                if line.find('*NSET, NSET=NODES_Z0') != -1:
                    nset_z0 = True
                if line.find('*STEP') != -1:
                    step = True
                if line.find('NODES_T, 3, 3, -0.002') != -1:
                    bc = True

        # read and check input CALCULIX file: read two specific file lines
        with open("foo.inp", 'r') as f:
            line_numbers = [99, 860]
            lines = []
            for i, line in enumerate(f):
                if i in line_numbers:
                    lines.append(line.strip())

        self.assertEqual(node, True)
        self.assertEqual(element, True)
        self.assertEqual(elset_z1, True)
        self.assertEqual(nset_z0, True)
        self.assertEqual(step, True)
        self.assertEqual(bc, True)

        # assert one node and one node_set line
        self.assertEqual(list(map(float, lines[0].split(","))), [359.0, 0.84, 1.2, 0.24])
        self.assertEqual(list(map(int, lines[1].strip(',').split(','))), [37, 64, 92, 126])

if __name__ == '__main__':
    unittest.main()