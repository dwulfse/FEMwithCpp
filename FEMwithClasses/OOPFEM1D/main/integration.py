import numpy as np
from scipy.integrate import dblquad

def integrate_over_triangle(f, vertices):
    """
    Integrate a function f(x,y) over a triangle defined by 3 vertices.
    
    Parameters:
      f: function f(x, y)
      vertices: list or array of three vertices, each [x, y]
      
    Returns:
      integral: the computed integral value,
      error: an estimate of the integration error.
    """
    v0 = np.array(vertices[0])
    v1 = np.array(vertices[1])
    v2 = np.array(vertices[2])
    
    # Compute the Jacobian (twice the area) of the affine mapping.
    jacobian = np.abs(np.linalg.det(np.column_stack((v1 - v0, v2 - v0))))
    
    def mapped_f(s, r):
        # Affine mapping: (r, s) in reference triangle (r in [0,1], s in [0,1-r])
        pt = v0 + r * (v1 - v0) + s * (v2 - v0)
        return f(pt[0], pt[1]) * jacobian

    # Integrate with r from 0 to 1 and s from 0 to (1 - r)
    integral, error = dblquad(mapped_f, 0, 1, lambda r: 0, lambda r: 1 - r)
    return integral, error

# Example usage:
if __name__ == '__main__':
    # Test with a constant function f(x,y)=1.
    # For a triangle with vertices (0,0), (1,0), (0,1), the area is 0.5.
    vertices = [[0.5, 0], [0, 0], [0.25, 0.25]]
    # f_const = lambda x, y: 1.0
    # integral, error = integrate_over_triangle(f_const, vertices)
    # print("Integral over triangle (expected 0.5):", integral)
    # print("Error estimate:", error)

    # You can also test with a non-constant function.
    f_test = lambda x, y: 2.0 * np.pi ** 2 * np.sin(np.pi * x) * np.sin(np.pi * y)
    integral2, error2 = integrate_over_triangle(f_test, vertices)
    print("\nIntegral of sin(pi*x)*sin(pi*y) over triangle:", integral2)
    print("Error estimate:", error2)
