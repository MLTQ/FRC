"""Regular-axis material flux labels and their geometric chain rule."""
import numpy as np
from scipy.sparse import coo_matrix


class MaterialFlux:
    def __init__(self, domain, initial_a):
        self.domain = domain
        self.labels = np.asarray(initial_a,float).copy()
        if self.labels.shape != (domain.fluid_node_count,) or not np.all(np.isfinite(self.labels)):
            raise ValueError("Finite initial a=psi/r² required, including its regular axis limit")
        axis_ids = np.flatnonzero(domain.axis)
        self.axis_ids = axis_ids
        rows,cols,values = [],[],[]
        for row,index in enumerate(axis_ids):
            if index == 0:
                adjacent = domain.ids[1,1:-1]
                weights = np.sin(domain.theta[1:-1])**2
                weights /= weights.sum()
            else:
                i,j = np.argwhere(domain.ids[:domain.nr+1] == index)[0]
                adjacent = np.array([domain.ids[i,1 if j == 0 else -2]])
                weights = np.ones(1)
            for node,weight in zip(adjacent,weights):
                rows.append(row);cols.append(node)
                values.append(weight/domain.fluid_reference[node,0])
        self.axis_stretch = coo_matrix((values,(rows,cols)),shape=(len(axis_ids),domain.fluid_node_count)).tocsr()

    def coefficients(self, coordinates):
        d = self.domain
        d.coordinates(coordinates)  # Validate before division by material radius.
        radius = coordinates[:,0]
        scale = np.ones(d.fluid_node_count)
        scale[~d.axis] = (d.fluid_reference[~d.axis,0]/radius[~d.axis])**2
        stretch = self.axis_stretch@radius
        if np.any(stretch <= 0) or not np.all(np.isfinite(stretch)):
            raise ValueError("Positive regular-axis radial stretch is required")
        scale[d.axis] = 1/stretch**2
        return self.labels*scale

    def pullback(self, coordinates, coefficient_load):
        """Radial derivative of load dot advected a, without changing labels."""
        d = self.domain;radius = coordinates[:,0]
        a = self.coefficients(coordinates)
        result = np.zeros_like(radius)
        result[~d.axis] = -2*a[~d.axis]*coefficient_load[~d.axis]/radius[~d.axis]
        stretch = self.axis_stretch@radius
        result += self.axis_stretch.T@(-2*a[d.axis]*coefficient_load[d.axis]/stretch)
        return result
