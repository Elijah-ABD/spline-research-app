**Catmull-Rom Spline Matrix**
```
| 0 -1/2 1 -1/2 |
| 1 0 -5/2 3/2 |
| 0 1/2 2 -3/2 |
| 0 0 1/2 1/2 |
```
This matrix is used to compute the control points for a Catmull-Rom spline, given the positions of the four surrounding points (p0, p1, p2, and p3).
**Centripetal Catmull-Rom Spline Matrix**
The Centripetal Catmull-Rom spline matrix is a variant of the Catmull-Rom spline matrix that produces a more uniform parameterization:
```
| 0 -1/2 1 -1/2 |
| 1 0 -3/2 5/4 |
| 0 1/2 1/2 -1/4 |
| 0 0 1/4 1/4 |
```
This matrix is used to compute the control points for a Centripetal Catmull-Rom spline, which is a variant of the Catmull-Rom spline that has a more uniform parameterization.
**Chordal Catmull-Rom Spline Matrix**
The Chordal Catmull-Rom spline matrix is another variant of the Catmull-Rom spline matrix that produces a more uniform parameterization:
```
| 0 -1/2 1 -1/2 |
| 1 0 -2 3/2 |
| 0 1/2 1 -1/2 |
| 0 0 1/2 1/2 |
```
This matrix is used to compute the control points for a Chordal Catmull-Rom spline, which is a variant of the Catmull-Rom spline that has a more uniform parameterization based on

q0 = -0.5 * tCu + tSq - 0.5 * t;
q1 = 1.5 * tCu - 2.5 * tSq + 1;
q2 = -1.5 * tCu + 2 * tSq + 0.5 * t;
q3 = 0.5 * tCu - 0.5 * tSq;
double tX = coords[p0].Item1 * q0 + coords[p1].Item1 * q1 + coords[p2].Item1 * q2 + coords[p3].Item1 * q3;
double tY = coords[p0].Item2 * q0 + coords[p1].Item2 * q1 + coords[p2].Item2 * q2 + coords[p3].Item2 * q3;
// ... (rest of the method remains the same)

// Centripetal coefficients
double alpha = 0.5; // adjust this value to change the shape of the curve
double t2 = t * t;
double t3 = t2 * t;
q0 = (1 - alpha) * (1 - t3) + alpha * (1 - 3 * t2 + 2 * t3);
q1 = (1 - alpha) * (3 * t - 2 * t2 - t3) + alpha * (4 * t - 4 * t2 + t3);
q2 = (1 - alpha) * (t2 - t3) + alpha * (t2 - t3);
q3 = (1 - alpha) * (t3) + alpha * (t3);
double tX = coords[p0].Item1 * q0 + coords[p1].Item1 * q1 + coords[p2].Item1 * q2 + coords[p3].Item1 * q3;
double tY = coords[p0].Item2 * q0 + coords[p1].Item2 * q1 + coords[p2].Item2 * q2 + coords[p3].Item2 * q3;