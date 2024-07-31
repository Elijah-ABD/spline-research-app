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