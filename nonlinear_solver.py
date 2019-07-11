import numpy as np


class NewtonRahpson:
    def __init__(self, nodes, elements):
        self._nodes = nodes
        self._elements = elements
        self.n_dof = len(nodes)*3

        self._u_bc = np.zeros(self.n_dof, dtype=bool)

        self.u = np.zeros(self.n_dof)
        self.f = np.zeros(self.n_dof)

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
            print iteration, self.convergence(du)
        print "The load step completed in", iteration, 'iterations'

    def convergence(self, du):
        if np.sum(du**2) == 0.:
            return 0
        return np.sqrt(np.dot(du, du))/(np.dot(self.u, self.u))

    def assembly(self):
        self._Kt = self._elements[0].Kt()
        self._fint = self._elements[0].f_int()

    # noinspection PyPep8Naming
    def _reduce(self, iteration):
        K_red = np.copy(self._Kt)
        bc_idx = np.where(self._u_bc)[0]
        R = self.f - self._fint
        for i in range(self.n_dof):
            for l in bc_idx:
                if iteration == 1:
                    du = self.u[l]
                else:
                    du = 0
                if i == l:
                    R[i] = du
                    K_red[l, l] = 1
                else:
                    R[i] -= K_red[i, l]*du
                    K_red[i, l] = 0
                    K_red[l, i] = 0
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
                self.u[idx] = value
