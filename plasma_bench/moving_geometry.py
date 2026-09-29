"""Material meridional mesh, massless magnetic mesh and curved matching sphere."""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from numpy.polynomial.legendre import Legendre


@dataclass
class Elements:
    field_ids: np.ndarray
    geometry_ids: np.ndarray
    field_shape: np.ndarray
    geometry_shape: np.ndarray
    field_gradient: np.ndarray
    geometry_gradient: np.ndarray
    points: np.ndarray
    weights: np.ndarray
    jacobian: np.ndarray

    def stiffness(self):
        return np.einsum("cq,cqna,cqma->cnm", self.weights*self.points[..., 0]**3,
                         self.field_gradient, self.field_gradient)

    def shape_gradient(self, values):
        """Derivative of integral r³ |grad a|² at fixed field coefficients."""
        gradient = np.einsum("cn,cqna->cqa", values[self.field_ids], self.field_gradient)
        square = np.sum(gradient*gradient, axis=-1)
        dot = np.einsum("cqa,cqna->cqn", gradient, self.geometry_gradient)
        radius = self.points[..., 0]
        result = radius[..., None, None]**3*(square[..., None, None]*self.geometry_gradient-
                                             2*dot[..., None]*gradient[..., None, :])
        result[..., 0] += 3*radius[..., None]**2*self.geometry_shape*square[..., None]
        return np.einsum("cq,cqna->cna", self.weights, result)


class MovingDomain:
    """Q1 material cells inside a deformable sphere-like boundary; no exterior gas."""
    def __init__(self, radial_cells=8, angular_cells=16, *, radius=1., matching_radius=2.,
                 vacuum_layers=None, modes=8, quadrature=5, radial_grading=1.):
        vacuum_layers = radial_cells if vacuum_layers is None else vacuum_layers
        for name, value, minimum in (("radial cells", radial_cells, 2), ("angular cells", angular_cells, 4),
                                     ("vacuum layers", vacuum_layers, 2), ("modes", modes, 1), ("quadrature", quadrature, 2)):
            if isinstance(value, bool) or not isinstance(value, int) or value < minimum:
                raise ValueError(f"Invalid {name}")
        if not np.isfinite(radius+matching_radius) or not 0 < radius < matching_radius:
            raise ValueError("Positive plasma radius must be inside the matching sphere")
        if modes > angular_cells:
            raise ValueError("Matching angular resolution must resolve the retained modes")
        if not np.isscalar(radial_grading) or not np.isfinite(radial_grading) or radial_grading <= 0:
            raise ValueError("Finite positive radial grading required")
        self.nr, self.nt, self.nv = radial_cells, angular_cells, vacuum_layers
        self.radius, self.matching_radius, self.modes = radius, matching_radius, modes
        self.quadrature, self.radial_grading = quadrature, float(radial_grading)
        theta = np.linspace(0., np.pi, angular_cells+1)
        shells = np.r_[0., np.linspace(radius/radial_cells, radius, radial_cells),
                       np.linspace(radius, matching_radius, vacuum_layers+1)[1:-1]]
        if radial_grading != 1.:
            material = radius*(1-(1-np.linspace(0.,1.,radial_cells+1))**radial_grading)
            exterior = radius+(matching_radius-radius)*np.linspace(0.,1.,vacuum_layers+1)[1:-1]**radial_grading
            shells = np.r_[material,exterior]
        if np.any(np.diff(np.r_[shells,matching_radius]) <= 0):
            raise ValueError("Radial grading collapses shell spacing at floating-point precision")
        self.theta, self.shells = theta, shells
        self.ids = np.zeros((len(shells), angular_cells+1), int)
        self.ids[1:] = np.arange(1, 1+(len(shells)-1)*(angular_cells+1)).reshape(len(shells)-1, angular_cells+1)
        points = np.stack((shells[1:, None]*np.sin(theta), -shells[1:, None]*np.cos(theta)), axis=-1)
        points[:, [0, -1], 0] = 0.
        self.reference = np.vstack((np.zeros((1, 2)), points.reshape(-1, 2)))
        self.fluid_node_count = 1+radial_cells*(angular_cells+1)
        self.node_count = len(self.reference)
        self.field_count = self.node_count+modes
        self.fluid_cell_count = radial_cells*angular_cells
        self.fluid_reference = self.reference[:self.fluid_node_count].copy()
        self.boundary = self.ids[radial_cells].copy()
        self.axis = self.fluid_reference[:, 0] == 0
        cells = []
        for i in range(len(shells)-1):
            for j in range(angular_cells):
                cells.append([self.ids[i,j], self.ids[i+1,j], self.ids[i+1,j+1], self.ids[i,j+1]])
        self.cells = np.asarray(cells)
        # Exterior geometric displacement depends only on the moving material boundary.
        rows, cols, values = list(range(self.fluid_node_count)), list(range(self.fluid_node_count)), [1.]*self.fluid_node_count
        for i in range(radial_cells+1, len(shells)):
            fraction = (matching_radius-shells[i])/(matching_radius-radius)
            for j in range(angular_cells+1):
                rows.append(self.ids[i,j]); cols.append(self.boundary[j]); values.append(fraction)
        self.motion = coo_matrix((values,(rows,cols)), shape=(self.node_count,self.fluid_node_count)).tocsr()
        x, w = np.polynomial.legendre.leggauss(quadrature)
        xi, tau = np.meshgrid((x+1)/2, (x+1)/2, indexing="ij")
        self.qweight = (.25*w[:,None]*w[None,:]).ravel()
        xi, tau = xi.ravel(), tau.ravel()
        self.shape = np.stack(((1-xi)*(1-tau), xi*(1-tau), xi*tau, (1-xi)*tau), axis=-1)
        self.derivative = np.stack((np.stack((-(1-tau), 1-tau, tau, -tau), axis=-1),
                                    np.stack((-(1-xi), -xi, xi, 1-xi), axis=-1)), axis=-1)
        self._outer_basis(xi, tau)
        ell = np.arange(1, modes+1)
        self.exterior_diagonal = (ell+2)*matching_radius**3
        self.physical_exterior_diagonal = ell*matching_radius**3

    def _outer_basis(self, xi, tau):
        width = np.pi/self.nt
        angle = self.theta[:-1,None]+width*tau[None,:]
        mu = -np.cos(angle)
        basis, derivative = [], []
        for ell in range(1, self.modes+1):
            polynomial = Legendre.basis(ell).deriv()
            norm = np.sqrt(2*ell*(ell+1)/(2*ell+1))
            basis.append(polynomial(mu)/norm)
            derivative.append(polynomial.deriv()(mu)*np.sin(angle)*width/norm)
        basis, derivative = np.stack(basis,axis=-1), np.stack(derivative,axis=-1)
        inner_shape = np.stack(((1-xi)*(1-tau),(1-xi)*tau),axis=-1)
        inner_derivative = np.stack((np.stack((-(1-tau),-tau),axis=-1),
                                     np.stack((-(1-xi),1-xi),axis=-1)),axis=-1)
        self.outer_geometry_shape = np.broadcast_to(inner_shape,(self.nt,*inner_shape.shape))
        self.outer_geometry_derivative = np.broadcast_to(inner_derivative,(self.nt,*inner_derivative.shape))
        self.outer_shape = np.concatenate((self.outer_geometry_shape,xi[None,:,None]*basis),axis=-1)
        modal_derivative = np.stack((basis,xi[None,:,None]*derivative),axis=-1)
        self.outer_derivative = np.concatenate((self.outer_geometry_derivative,modal_derivative),axis=-2)
        arc = self.matching_radius*np.stack((np.sin(angle),-np.cos(angle)),axis=-1)
        arc_derivative = self.matching_radius*width*np.stack((np.cos(angle),np.sin(angle)),axis=-1)
        self.outer_fixed = xi[None,:,None]*arc
        self.outer_fixed_derivative = np.stack((arc,xi[None,:,None]*arc_derivative),axis=-1)
        inner = np.column_stack((self.ids[-1,:-1],self.ids[-1,1:]))
        self.outer_geometry_ids = inner
        self.outer_field_ids = np.column_stack((inner,np.broadcast_to(np.arange(self.node_count,self.field_count),(self.nt,self.modes))))

    def coordinates(self, fluid):
        fluid = np.asarray(fluid,float)
        if fluid.shape != self.fluid_reference.shape or not np.all(np.isfinite(fluid)):
            raise ValueError("Finite material-node coordinates must match the domain")
        if np.any(fluid[self.axis,0] != 0) or np.any(fluid[~self.axis,0] <= 0):
            raise ValueError("Axis must stay at r=0 and other material nodes at r>0")
        if np.any(np.linalg.norm(fluid,axis=-1) >= self.matching_radius):
            raise ValueError("The fixed matching sphere must enclose every material node")
        return self.reference+self.motion@(fluid-self.fluid_reference)

    def elements(self, fluid):
        coordinates = self.coordinates(fluid)
        corner_derivative = np.array([[[-1.,-1.],[1.,0.],[0.,0.],[0.,1.]],
                                      [[-1.,0.],[1.,-1.],[0.,1.],[0.,0.]],
                                      [[0.,0.],[0.,-1.],[1.,1.],[-1.,0.]],
                                      [[0.,-1.],[0.,0.],[1.,0.],[-1.,1.]]])
        corner_jac = np.einsum("qnb,cna->cqab",corner_derivative,coordinates[self.cells])
        corner_det = np.linalg.det(corner_jac)
        allowed_zero = np.zeros_like(corner_det,dtype=bool)
        allowed_zero[:self.nt,[0,3]] = True  # Intentionally collapsed center edge.
        if np.any(corner_det < 0) or np.any((corner_det == 0)&~allowed_zero):
            raise ValueError("Nonpositive cell-corner Jacobian; material mesh cannot be repaired silently")
        result = []
        for ids, geometry_ids, shape, deriv, geo_shape, geo_deriv, fixed, fixed_deriv in (
            (self.cells,self.cells,self.shape,self.derivative,self.shape,self.derivative,0.,0.),
            (self.outer_field_ids,self.outer_geometry_ids,self.outer_shape,self.outer_derivative,
             self.outer_geometry_shape,self.outer_geometry_derivative,self.outer_fixed,self.outer_fixed_derivative)):
            count = len(ids)
            shape = np.broadcast_to(shape,(count,*shape.shape[-2:]))
            deriv = np.broadcast_to(deriv,(count,*deriv.shape[-3:]))
            geo_shape = np.broadcast_to(geo_shape,(count,*geo_shape.shape[-2:]))
            geo_deriv = np.broadcast_to(geo_deriv,(count,*geo_deriv.shape[-3:]))
            points = np.einsum("cqn,cna->cqa",geo_shape,coordinates[geometry_ids])+fixed
            jac = np.einsum("cqnb,cna->cqab",geo_deriv,coordinates[geometry_ids])+fixed_deriv
            det = np.linalg.det(jac)
            if np.any(det <= 0) or np.any(points[...,0] <= 0) or not np.all(np.isfinite(det)):
                raise ValueError("Inverted or nonphysical material/vacuum mesh; no repair applied")
            inverse = np.linalg.inv(jac)
            gradient = np.einsum("cqnb,cqba->cqna",deriv,inverse)
            geo_gradient = np.einsum("cqnb,cqba->cqna",geo_deriv,inverse)
            result.append(Elements(ids,geometry_ids,shape,geo_shape,gradient,geo_gradient,
                                   points,self.qweight[None,:]*det,det))
        return tuple(result)

    def stiffness(self, elements):
        rows,cols,values = [],[],[]
        for batch in elements:
            ids = batch.field_ids; width = ids.shape[1]
            rows.append(np.repeat(ids,width,axis=1).ravel())
            cols.append(np.tile(ids,(1,width)).ravel())
            values.append(batch.stiffness().ravel())
        modal = np.arange(self.node_count,self.field_count)
        rows.append(modal);cols.append(modal);values.append(self.exterior_diagonal)
        return coo_matrix((np.concatenate(values),(np.concatenate(rows),np.concatenate(cols))),
                          shape=(self.field_count,self.field_count)).tocsc()

    def volumes(self, elements):
        batch = elements[0]; n = self.fluid_cell_count
        return 2*np.pi*np.sum(batch.weights[:n]*batch.points[:n,:,0],axis=1)

    def volume_gradients(self, elements):
        batch = elements[0]; n = self.fluid_cell_count
        gradient = batch.points[:n,:,0,None,None]*batch.geometry_gradient[:n]
        gradient = gradient.copy(); gradient[...,0] += batch.geometry_shape[:n]
        return 2*np.pi*np.einsum("cq,cqna->cna",batch.weights[:n],gradient)
