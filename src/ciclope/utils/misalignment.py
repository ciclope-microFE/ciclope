import matplotlib.pyplot as plt

import numpy as np
from numpy.linalg import qr, svd, eig, norm, solve

import warnings

from skimage.measure import find_contours
from scipy.spatial import ConvexHull
from scipy.optimize import minimize_scalar, least_squares

#-------------------------------------------------------------------------------------

def point_plane_signed_distance(P, origin, normal):
    """
    Compute the signed distance from a point to a plane.

    The signed distance is positive if the point lies in the 
    direction of the plane normal, negative otherwise.

    Parameters
    ----------
    P : array_like
        Coordinates of the point (3-element array).
    origin : array_like
        A point on the plane (3-element array).
    normal : array_like
        The normal vector of the plane (3-element array).

    Returns
    -------
    float
        The signed distance from the point to the plane.
    """
    return np.dot(P - origin, normal)

#-------------------------------------------------------------------------------------

def best_fit_segment_3d(points):
    """
    Compute the best-fit line segment through a set of 3D points.

    This function fits a line to the input points in a least-squares sense 
    using Singular Value Decomposition (SVD). It returns the endpoints of 
    the segment that spans the projection of the input points onto the 
    best-fit line, as well as the unit direction vector.

    Parameters
    ----------
    points : array_like, shape (N, 3)
        Array of 3D points to fit. Each row represents a point.

    Returns
    -------
    p1 : ndarray, shape (3,)
        First endpoint of the best-fit segment.
    p2 : ndarray, shape (3,)
        Second endpoint of the best-fit segment.
    direction : ndarray, shape (3,)
        Unit direction vector of the best-fit line.
    """
    
    points = np.asarray(points)
    centroid = np.mean(points, axis=0)

    # Principal direction with SVD
    _, _, vh = np.linalg.svd(points - centroid)
    direction = vh[0]  # unit direction vector

    # Project points onto the line
    diffs = points - centroid
    t_vals = np.dot(diffs, direction)

    # Find min and max projections
    t_min, t_max = t_vals.min(), t_vals.max()
    p1 = centroid + t_min * direction
    p2 = centroid + t_max * direction

    return p1, p2, direction

#-------------------------------------------------------------------------------------

def signed_point_to_ellipse_distance(P, C, a, b, alpha):
    """
    Compute the signed minimum distance between a point and a rotated ellipse in 2D.

    The distance is positive if the point lies outside the ellipse, negative if inside.
    The ellipse can be arbitrarily rotated about its center.

    Parameters
    ----------
    P : array_like, shape (2,)
        Coordinates of the point in the XY plane.
    C : array_like, shape (2,)
        Center of the ellipse.
    a : float
        Semi-major axis length.
    b : float
        Semi-minor axis length.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.

    Returns
    -------
    float
        Signed Euclidean distance from the point to the ellipse.
        Positive if the point is outside, negative if inside.
    """
    
    # Convert P and C to numpy array
    P = np.asarray(P)
    C = np.asarray(C)

    # Rotate and translate the point in the space of the unrotated ellipse
    cos_a, sin_a = np.cos(-alpha), np.sin(-alpha)
    R = np.array([[cos_a, -sin_a],
                  [sin_a, cos_a]])
    
    P_rot = R @ (P - C)

    # Square distance function from a point on a parametric ellipse
    def distance_squared(theta):
        x = a * np.cos(theta)
        y = b * np.sin(theta)
        return (x - P_rot[0])**2 + (y - P_rot[1])**2

    # Minimize the distance squared to find the nearest point on the ellipse
    res = minimize_scalar(distance_squared, bounds=(0, 2*np.pi), method='bounded')

    # Calculate the minimum point on the ellipse (not rotated)
    theta_min = res.x
    x_min = a * np.cos(theta_min)
    y_min = b * np.sin(theta_min)

    # Rotate point found in original space
    cos_a, sin_a = np.cos(alpha), np.sin(alpha)
    R_inv = np.array([[cos_a, -sin_a],
                      [sin_a, cos_a]])
    closest_point = R_inv @ np.array([x_min, y_min]) + C

    # Euclidean distance from original point
    euclideanDist=np.sqrt((P[0] - closest_point[0])**2 + (P[1] - closest_point[1])**2)

    centerClosDist=np.sqrt((C[0] - closest_point[0])**2 + (C[1] - closest_point[1])**2)
    centerPDist=np.sqrt((C[0] - P[0])**2 + (C[1] - P[1])**2)

    if centerClosDist < centerPDist :
        return euclideanDist
    else :
        return -euclideanDist
    
#-------------------------------------------------------------------------------------

def fit_ellipse(x, linear=False, constraint='bookstein', maxits=200, tol=1e-5):
    """
    Fit an ellipse to a set of 2D points.

    This function fits an ellipse in least-squares sense using either a 
    linear algebraic method (Bookstein or Trace constraint) and optionally 
    refines the result with a nonlinear Gauss-Newton optimization.

    Parameters
    ----------
    x : array_like, shape (N, 2) or (2, N)
        2D coordinates of the points to fit. Can be provided as N×2 or 2×N array.
    linear : bool, optional
        If True, only the linear fit is performed. If False (default), a nonlinear 
        refinement is applied after the linear estimate.
    constraint : {'bookstein', 'trace'}, optional
        The type of constraint for the linear fit. 'bookstein' uses the Bookstein 
        constraint, 'trace' uses the Trace constraint.
    maxits : int, optional
        Maximum number of Gauss-Newton iterations for nonlinear refinement.
    tol : float, optional
        Convergence tolerance for nonlinear refinement.

    Returns
    -------
    z : ndarray, shape (2,)
        Coordinates of the ellipse center.
    a : float
        Semi-major axis length.
    b : float
        Semi-minor axis length.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.
    """

    # Prepare input
    x = np.asarray(x)
    if x.ndim != 2 or x.shape[0] not in (2, x.shape[1]):
        if x.shape[1] == 2:
            x = x.T
        else:
            raise ValueError("Input must be 2xN or Nx2 array")
    if x.shape[1] < 6:
        raise ValueError("At least 6 points required to compute fit")

    # Remove centroid for numerical stability
    centroid = x.mean(axis=1, keepdims=True)
    x_cent = x - centroid

    # Linear fit
    if constraint == 'bookstein':
        z, a, b, alpha = _fit_bookstein(x_cent)
    elif constraint == 'trace':
        z, a, b, alpha = _fit_ggk(x_cent)
    else:
        raise ValueError("Invalid constraint: choose 'bookstein' or 'trace'")

    # Nonlinear refinement
    if not linear:
        try:
            z_ref, a_ref, b_ref, alpha_ref, converged = _fit_nonlinear(
                x_cent, z, a, b, alpha, maxits, tol)
            if converged:
                z, a, b, alpha = z_ref, a_ref, b_ref, alpha_ref
            #else:
                #warnings.warn("Gauss-Newton did not converge, using linear estimate")
        except Exception as e:
            warnings.warn(f"Nonlinear fit failed: {e}")

    # Add back centroid
    z = z.flatten() + centroid.flatten()
    return z, a, b, alpha

#-------------------------------------------------------------------------------------

def fit_ellipse_filtered(pts):
    """
    Fit an ellipse (or circle) to 2D points with outlier filtering.

    This function first fits an ellipse to the input points, then iteratively filters 
    out points that deviate beyond a threshold percentage of the average semi-axis length.
    The fit is repeated up to 4 times with progressively tighter thresholds.

    Parameters
    ----------
    pts : array_like, shape (N, 2)
        Array of 2D points to fit, as [[x1, y1], [x2, y2], ...].

    Returns
    -------
    C : ndarray, shape (2,)
        Coordinates of the ellipse center after filtering.
    a : float
        Semi-major axis length after filtering.
    b : float
        Semi-minor axis length after filtering.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.
    """
    
    threshold_percent = 5
    # First fit of the circle
    C, a, b, alpha = fit_ellipse(pts)
    
    for i in range(0,4):
        # Calculate distances of points from the centre
        distances=np.array([signed_point_to_ellipse_distance(p, C, a, b, alpha) for p in pts])
        
        # Distance threshold (5% of mean semi-axis)
        medium_semi_ax=(a+b)/2.0 
        

        threshold = medium_semi_ax* ((threshold_percent-i) / 100)
        
        # Filter the points that are within the threshold
        filtered_pts = pts[distances > -threshold]

        # Second fit of the circle with filtered points
        C, a, b, alpha = fit_ellipse(filtered_pts)
        
    return C, a, b, alpha

#-------------------------------------------------------------------------------------

def _fit_bookstein(x):
    """
    Fit an ellipse to 2D points using the Bookstein constraint.

    This is a helper function that performs a direct least-squares fit of 
    an ellipse using the Bookstein constraint, which enforces that the 
    sum of squared ellipse coefficients equals one for numerical stability.

    Parameters
    ----------
    x : ndarray, shape (2, N)
        2D points as a 2×N array.

    Returns
    -------
    z : ndarray, shape (2,)
        Coordinates of the ellipse center.
    a : float
        Semi-major axis length.
    b : float
        Semi-minor axis length.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.

    Notes
    -----
    This function is typically called internally by `fit_ellipse` and is not 
    intended to be used directly.
    """

    m = x.shape[1]
    x1 = x[0, :]
    x2 = x[1, :]
    B = np.column_stack([x1, x2, np.ones(m), x1**2,
                         np.sqrt(2)*x1*x2, x2**2])
    Q, R = qr(B, mode='reduced')
    R11 = R[:3, :3]
    R12 = R[:3, 3:6]
    R22 = R[3:6, 3:6]
    U, S, Vt = svd(R22)
    w = Vt.T[:, 2]
    v = -solve(R11, R12.dot(w))
    # Build quadratic form
    A = np.array([[w[0], w[1]/np.sqrt(2)], [w[1]/np.sqrt(2), w[2]]])
    bv = v[:2]
    c = v[2]
    return _conic2parametric(A, bv, c)

#-------------------------------------------------------------------------------------

def _fit_ggk(x):
    """
    Fit an ellipse to 2D points using the Trace constraint (GGK method).

    This helper function performs a direct least-squares fit of an ellipse 
    using the Trace constraint, which enforces that the sum of the ellipse's 
    quadratic form eigenvalues equals one for numerical stability.

    Parameters
    ----------
    x : ndarray, shape (2, N)
        2D points as a 2×N array.

    Returns
    -------
    z : ndarray, shape (2,)
        Coordinates of the ellipse center.
    a : float
        Semi-major axis length.
    b : float
        Semi-minor axis length.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.

    Notes
    -----
    This function is used internally by `fit_ellipse` when the 'trace' 
    constraint is selected and is not intended to be used directly.
    """

    m = x.shape[1]
    x1 = x[0, :]
    x2 = x[1, :]
    B = np.column_stack([2*x1*x2, x2**2 - x1**2, x1, x2, np.ones(m)])
    v = solve(B, -x1**2)
    A = np.array([[1 - v[1], v[0]], [v[0], v[1]]])
    bv = v[2:4]
    c = v[4]
    return _conic2parametric(A, bv, c)

#-------------------------------------------------------------------------------------

def _fit_nonlinear(x, z0, a0, b0, alpha0, maxits, tol):
    """
    Refine ellipse parameters by nonlinear least squares optimization.

    This function performs Gauss-Newton iterations to minimize the geometric
    error between the ellipse model and the given 2D points, refining initial
    estimates of the ellipse center, axes, and rotation.

    Parameters
    ----------
    x : ndarray, shape (2, N)
        2D points as a 2×N array.
    z0 : ndarray, shape (2,)
        Initial estimate of the ellipse center.
    a0 : float
        Initial estimate of the semi-major axis length.
    b0 : float
        Initial estimate of the semi-minor axis length.
    alpha0 : float
        Initial estimate of the ellipse rotation angle in radians.
    maxits : int
        Maximum number of Gauss-Newton iterations.
    tol : float
        Convergence tolerance based on relative update size.

    Returns
    -------
    z : ndarray, shape (2,)
        Refined ellipse center coordinates.
    a : float
        Refined semi-major axis length.
    b : float
        Refined semi-minor axis length.
    alpha : float
        Refined rotation angle in radians.
    converged : bool
        True if the optimization converged within the iteration limit.

    Notes
    -----
    This function is typically used internally by `fit_ellipse` during 
    nonlinear refinement and is not intended for direct use.
    """

    m = x.shape[1]
    # initial phase phi0
    Q0 = np.array([[np.cos(alpha0), -np.sin(alpha0)],
                   [np.sin(alpha0), np.cos(alpha0)]])
    diffs = x - z0.reshape(2,1)
    vals = Q0.T.dot(diffs)
    phi0 = np.arctan2(vals[1,:], vals[0,:])
    # state vector u = [phi..., alpha, a, b, z1, z2]
    u = np.concatenate([phi0, [alpha0, a0, b0], z0.flatten()])
    converged = False
    for _ in range(maxits):
        f, J = _sys_fn(u, x)
        try:
            h = solve(J, -f)
        except np.linalg.LinAlgError:
            break
        u += h
        if norm(h, np.inf)/norm(u, np.inf) < tol:
            converged = True
            break
    alpha = u[m]
    a = u[m+1]
    b = u[m+2]
    z = u[m+3:m+5]
    return z, a, b, alpha, converged

#-------------------------------------------------------------------------------------

def _sys_fn(u, x):
    """
    Compute residuals and Jacobian matrix for ellipse fitting nonlinear system.

    Given the current parameter estimates, this function calculates the residuals
    between the observed points and their projections on the ellipse, along with 
    the Jacobian matrix of partial derivatives needed for Gauss-Newton optimization.

    Parameters
    ----------
    u : ndarray, shape (m + 5,)
        State vector containing:
        - phi (m,): parameter angles for each point on the ellipse
        - alpha: rotation angle of the ellipse
        - a: semi-major axis length
        - b: semi-minor axis length
        - z (2,): ellipse center coordinates
    x : ndarray, shape (2, m)
        Observed 2D points as a 2×m array.

    Returns
    -------
    f : ndarray, shape (2*m,)
        Residual vector concatenating x - ellipse model points for all data.
    J : ndarray, shape (2*m, m + 5)
        Jacobian matrix of partial derivatives of residuals with respect to u.

    Notes
    -----
    This function is used internally by `_fit_nonlinear` during Gauss-Newton iterations
    and is not intended for external use.
    """

    m = x.shape[1]
    phi = u[:m]
    alpha = u[m]
    a = u[m+1]
    b = u[m+2]
    z = u[m+3:m+5]
    ca, sa = np.cos(alpha), np.sin(alpha)
    Q = np.array([[ca, -sa],[sa, ca]])
    Qdot = np.array([[-sa, -ca],[ca, -sa]])
    f = np.zeros(2*m)
    J = np.zeros((2*m, m+5))
    for i in range(m):
        idx = slice(2*i, 2*i+2)
        cphi, sphi = np.cos(phi[i]), np.sin(phi[i])
        f[idx] = x[:, i] - z - Q.dot([a*cphi, b*sphi])
        # Jacobian
        J[idx, i] = -Q.dot([-a*sphi, b*cphi])
        J[idx, m] = -Qdot.dot([a*cphi, b*sphi])
        J[idx, m+1] = -Q.dot([cphi, 0])
        J[idx, m+2] = -Q.dot([0, sphi])
        J[idx, m+3:m+5] = -np.eye(2)
    return f, J

#-------------------------------------------------------------------------------------

def _conic2parametric(A, bv, c):
    """
    Convert conic parameters of an ellipse to parametric form.

    Given the quadratic form parameters of an ellipse, this function computes
    the center coordinates, semi-major and semi-minor axes lengths, and rotation angle.

    Parameters
    ----------
    A : ndarray, shape (2, 2)
        Symmetric matrix representing the quadratic form of the ellipse.
    bv : ndarray, shape (2,)
        Linear terms vector of the ellipse equation.
    c : float
        Constant term of the ellipse equation.

    Returns
    -------
    t : ndarray, shape (2, 1)
        Coordinates of the ellipse center.
    a : float
        Semi-major axis length.
    b : float
        Semi-minor axis length.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.

    Raises
    ------
    ValueError
        If the conic does not represent a valid ellipse (e.g., eigenvalues product ≤ 0).

    Notes
    -----
    The input conic equation is assumed to be of the form:
    x.T @ A @ x + bv.T @ x + c = 0
    """

    # Diagonalize A: A = Q.T D Q
    eigvals, eigvecs = eig(A)
    if np.prod(eigvals) <= 0:
        raise ValueError("Linear fit did not produce an ellipse")
    Q = eigvecs.T
    D = np.diag(eigvals)
    t = -0.5 * solve(A, bv)
    c_h = t.T.dot(A).dot(t) + bv.T.dot(t) + c
    a = np.sqrt(-c_h / D[0,0])
    b = np.sqrt(-c_h / D[1,1])
    alpha = np.arctan2(Q[0,1], Q[0,0])
    return t.reshape(2,1), a, b, alpha

#--------------------------------------------------------------------

def plot_ellipse_3d(C, a, b, alpha, z_level, *, ax=None, line_spec=None, npts=100):
    """
    Plot a 2D ellipse in 3D space at a specified z-level.

    Parameters
    ----------
    C : array_like, shape (2,)
        Center coordinates of the ellipse in the XY plane.
    a : float
        Semi-major axis length of the ellipse.
    b : float
        Semi-minor axis length of the ellipse.
    alpha : float
        Rotation angle of the ellipse in radians, counterclockwise from the X-axis.
    z_level : float
        The Z coordinate at which to plot the ellipse.
    ax : matplotlib.axes._subplots.Axes3DSubplot, optional
        Existing 3D axes to plot on. If None, a new 3D axes is created.
    line_spec : str, optional
        Line specification string for plotting style (e.g., 'r-', 'b--').
    npts : int, optional
        Number of points to use for plotting the ellipse (default is 100).

    Returns
    -------
    line : matplotlib.lines.Line3D
        The matplotlib 3D line object representing the ellipse.

    Notes
    -----
    The ellipse is parameterized and rotated in the XY plane, then plotted
    at a constant Z height in 3D space.
    """

    # Obtain or create 3D axes
    if ax is None:
        fig = plt.gcf()
        ax = fig.add_subplot(111, projection='3d')

    # Parameter vector
    t = np.linspace(0, 2 * np.pi, npts)

    # Rotation matrix
    ca, sa = np.cos(alpha), np.sin(alpha)
    Q = np.array([[ca, -sa], [sa, ca]])

    # Ellipse points in 2D
    points_2d = Q.dot(np.vstack((a * np.cos(t), b * np.sin(t)))) + np.reshape(C, (2, 1))
    xs = points_2d[0, :]
    ys = points_2d[1, :]
    zs = np.full_like(xs, z_level)

    #print(f"Plotting ellipse at z={z_level} with center={C}, a={a}, b={b}, alpha={alpha}")
    #print(f"Ellipse X range: {xs.min()} - {xs.max()}, Y range: {ys.min()} - {ys.max()}")

    # Plot
    if line_spec:
        line, = ax.plot(xs, ys, zs, line_spec)
    else:
        line, = ax.plot(xs, ys, zs)

    return line

#-------------------------------------------------------------------------------------------

def ellipse_slice_uCT_trab_bone(mask2d, z_level, bin_thresh=117, ax=None):
    """
    Process a microCT image slice, fit an ellipse to the trabecular bone cross-section,
    and plot it in 3D at height z_level.
    
    This function extracts contours from a 2D binary mask slice, computes the 
    convex hull of the contour points, fits an ellipse to the hull points using 
    filtered fitting, and plots the resulting ellipse at the given z-level.

    Parameters
    ----------
    mask2d : ndarray, shape (H, W)
        2D binary mask slice where trabecular bone is represented.
    z_level : float
        The z-coordinate level for plotting the ellipse in 3D.
    bin_thresh : int, optional
        Threshold value for binarization (default is 117, currently unused).
    ax : matplotlib.axes._subplots.Axes3DSubplot, optional
        3D axes object to plot on. If None, plotting is skipped or uses current axes.

    Returns
    -------
    C : ndarray, shape (2,)
        Center coordinates of the fitted ellipse.
    a : float
        Semi-major axis length of the fitted ellipse.
    b : float
        Semi-minor axis length of the fitted ellipse.
    hull : scipy.spatial.ConvexHull
        Convex hull object computed from contour points.
    pts : ndarray, shape (N, 2)
        Points used for ellipse fitting (convex hull vertices).

    Notes
    -----
    If no contours are found in the mask slice, the function returns a dummy ellipse
    and prints a warning message.
    """
    
    # Find all contours in binary mask
    contours = find_contours(mask2d.astype(float), 0.5)
    if not contours:
        print(f"No contours found in sclice: {z_level}")
        return np.array([0, 0]), 1, 1  # dummy ellipse

    # Aggregate contour points
    coords = np.vstack(contours)
    # coords: rows (y), cols (x)

    # Compute convex hull of contour points
    hull = ConvexHull(coords)
    hull_pts = coords[hull.vertices]

    # Prepare points for ellipse fitting: x, y with origin at bottom-left
    ys = hull_pts[:, 0]
    xs = hull_pts[:, 1]
    pts = np.column_stack([xs, ys])
    
    # Fit ellipse
    C, a, b, alpha = fit_ellipse_filtered(pts)
    #print(f"Fitted ellipse: center={C}, a={a}, b={b}, alpha={alpha}")

    # Plot
    plot_ellipse_3d(C, a, b, alpha, z_level, ax=ax, line_spec='r--')

    return C, a, b, hull, pts

#-------------------------------------------------------------------------------------------

def calc_radius(center, points):
    """
    Calculate Euclidean distances from a center point to a set of 2D points.

    Parameters
    ----------
    center : array_like, shape (2,)
        Coordinates of the center point (x, y).
    points : ndarray, shape (N, 2)
        Array of 2D points with each row representing a point (x, y).

    Returns
    -------
    distances : ndarray, shape (N,)
        Array of Euclidean distances from the center to each point.
    """

    # distance from centre to each point
    return np.sqrt((points[:,0] - center[0])**2 + (points[:,1] - center[1])**2)

#-------------------------------------------------------------------------------------------

def residuals(params, points):
    """
    Compute residuals between distances from a center and a target radius.

    Parameters
    ----------
    params : array_like, shape (3,)
        Parameters array containing:
        - a : float, x-coordinate of the center
        - b : float, y-coordinate of the center
        - r : float, target radius
    points : ndarray, shape (N, 2)
        Array of 2D points to compare against the circle defined by the center and radius.

    Returns
    -------
    residuals : ndarray, shape (N,)
        Differences between the distance of each point from the center and the radius.
    """

    a, b, r = params
    return calc_radius([a,b], points) - r

#-------------------------------------------------------------------------------------------

def fit_circle(points):
    """
    Fit a circle to a set of 2D points using least squares optimization.

    Parameters
    ----------
    points : ndarray, shape (N, 2)
        Array of 2D points to fit the circle to.

    Returns
    -------
    a : float
        x-coordinate of the fitted circle's center.
    b : float
        y-coordinate of the fitted circle's center.
    r : float
        Radius of the fitted circle.

    Notes
    -----
    The initial guess for the circle center is the mean of the points,
    and the initial radius is the mean distance of the points from this center.
    The optimization minimizes the residuals between the distances of points
    from the center and the radius.
    """

    # Initial centre estimate as average of points
    center_estimate = np.mean(points, axis=0)
    # Initial radius estimate as the average of distances from the estimated centre
    r_estimate = np.mean(calc_radius(center_estimate, points))
    
    initial_params = [center_estimate[0], center_estimate[1], r_estimate]
    
    result = least_squares(residuals, initial_params, args=(points,))
    
    a, b, r = result.x
    return a, b, r

#-------------------------------------------------------------------------------------------

def fit_circle_filtered(pts, threshold_percent=5):
    """
    Fit a circle to 2D points with iterative outlier filtering based on distance thresholds.

    Parameters
    ----------
    pts : ndarray, shape (N, 2)
        Array of 2D points to fit the circle to.
    threshold_percent : float, optional
        Initial percentage of the radius used as a threshold for filtering points
        away from the fitted circle (default is 5%).

    Returns
    -------
    x : float
        x-coordinate of the fitted circle's center after filtering.
    y : float
        y-coordinate of the fitted circle's center after filtering.
    radius : float
        Radius of the fitted circle after filtering.

    Notes
    -----
    The fitting is performed iteratively:
    - Initially, the circle is fit to all points.
    - Points farther than a decreasing threshold percentage of the radius
      are filtered out.
    - The circle is refit to the filtered points for up to 4 iterations.
    """

    # First fit of the circle
    x, y, radius = fit_circle(pts)
    
    for i in range(0,4):
        # Calculate distances of points from the centre
        distances = np.sqrt((pts[:,0] - x)**2 + (pts[:,1] - y)**2)
        
        # Distance threshold (5% of radius)
        threshold = radius * ((threshold_percent-i) / 100)
        
        # Filter the points that are within the threshold
        filtered_pts = pts[distances >= (radius - threshold)]

        # Second fit of the circle with filtered points
        x, y, radius = fit_circle(filtered_pts)
        
    return x, y, radius

#-------------------------------------------------------------------------------------------

def circle_slice_uCT_trab_bone(mask2d, z_level, bin_thresh=1, ax=None):
    """
    Process a microCT image slice, fit an ellipse to the trabecular bone cross-section,
    and plot it in 3D at height z_level.

    Parameters
    ----------
    mask2d : ndarray, shape (M, N)
        2D binary mask representing a slice of trabecular bone.
    z_level : float or int
        The z-coordinate (height) at which the slice is located.
    bin_thresh : int, optional
        Threshold value used for binarization (default is 117).
    ax : matplotlib.axes._subplots.Axes3DSubplot, optional
        Matplotlib 3D axis for plotting the fitted circle.

    Returns
    -------
    C : list of float
        Center coordinates [x, y] of the fitted circle.
    r : float
        Radius of the fitted circle.
    r : float
        Radius repeated for compatibility with ellipse interface.
    hull : scipy.spatial.ConvexHull
        Convex hull object of the contour points.
    pts : ndarray, shape (K, 2)
        Array of points used for fitting the circle (convex hull vertices).

    Notes
    -----
    - Contours are extracted from the binary mask and used to compute the convex hull.
    - The circle is fit to the hull points with iterative outlier filtering.
    - The fitted circle is plotted on the specified z-level in 3D.
    """

    # Find all contours in binary mask
    contours = find_contours(mask2d.astype(float), 0.5)
    if not contours:
        print(f"No contours found in sclice: {z_level}")
        return np.array([0, 0]), 1, 1  # dummy ellipse

    # Aggregate contour points
    coords = np.vstack(contours)

    # Compute convex hull of contour points
    hull = ConvexHull(coords)
    hull_pts = coords[hull.vertices]

    # Prepare points for ellipse fitting: x, y with origin at bottom-left
    ys = hull_pts[:, 0]
    xs = hull_pts[:, 1]
    pts = np.column_stack([xs, ys])

    # Fit ellipse
    x, y, r = fit_circle_filtered(pts)
    #print(f"Fitted ellipse: center={C}, a={a}, b={b}, alpha={alpha}")
    C=[x,y]

    # Plot
    plot_ellipse_3d(C, r, r, 0, z_level, ax=ax, line_spec='r--')

    return C, r, r, hull, pts