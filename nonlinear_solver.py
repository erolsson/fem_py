import numpy as np


class NewtonRahpson:
    def __init__(self, nodes, elements):
        self._nodes = nodes
        self._elements = elements
        self.n_dof = len(nodes)*3
        self._u_bc = np.zeros(self.n_dof, dtype=bool)
        self._du = np.zeros(self.n_dof)

        self.u = np.zeros(self.n_dof)
        self._f = np.zeros(self.n_dof)

        self._Kt = np.zeros((self.n_dof, self.n_dof))
        self._fint = np.zeros(self.n_dof)

        self.tol = 1e-10
        self.assembly()

    # noinspection PyPep8Naming
    def solve(self):
        convergent = False
        iteration = 0
        while not convergent:
            K_red, R = self._reduce(iteration)
            du = np.linalg.solve(K_red, R)
            self.u += du
            self.update_elements()
            convergent = self.convergence(du) < self.tol
            iteration += 1
            print "iteration", iteration, ", residual", self.convergence(du)
        print "The load step completed in", iteration, 'iterations', self._elements[0]._gauss_points[0].material._fM

    @staticmethod
    def convergence(du):
        if np.sum(du**2) == 0.:
            return 0
        return np.sqrt(np.dot(du, du))

    def assembly(self):
        self._Kt = self._elements[0].Kt()
        self._fint = self._elements[0].f_int()

    # noinspection PyPep8Naming
    def _reduce(self, iteration):
        K_red = np.copy(self._Kt)
        R = self._f - self._fint
        du = np.copy(self._du)
        if iteration > 0:
            du *= 0
        for i in range(self.n_dof):
            rhs = 0
            for l in range(self.n_dof):
                if self._u_bc[l]:
                    rhs += K_red[i, l]*du[l]
                    if i == l:
                        K_red[i, l] = 1
                    else:
                        K_red[i, l] = 0
                        K_red[l, i] = 0
            if self._u_bc[i]:
                R[i] = du[i]
            else:
                R[i] -= rhs

        return K_red, R

    def update_elements(self):
        for e in self._elements:
            e.update(self.u)
        self._fint = self._elements[0].f_int()

    def set_displacement_bc(self, nodes, components, value=0):
        for node in nodes:
            for comp in components:
                idx = node*3 + comp
                self._u_bc[idx] = True
                self._du[idx] = value

    def set_force_bc(self, nodes, components, value=0):
        for node in nodes:
            for comp in components:
                idx = node*3 + comp
                self._f[idx] = value

    def reaction_forces(self):
        return self._fint
