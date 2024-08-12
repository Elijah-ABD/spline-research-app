using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Navigation;
using System.Windows.Shapes;
using MathNet.Numerics;

namespace WpfApp1
{
    public partial class MainWindow : System.Windows.Window
    {
        private List<Point> coords = new List<Point>();
        private CatmullRom cr = new CatmullRom();
        private float alpha = 0.5f;
        
        public MainWindow()
        {
            InitializeComponent();
        }

        // Event Handlers

        private void MouseClickHandler(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            coords.Add(point);
            Draw(point, 6, Brushes.LightGray);
                  
            if (coords.Count > 1)
            {
                SplineDrawing(Brushes.Orange);
            }
        }
public void DrawCone(Point p1, Point p2, double angleDegrees, Canvas canvas)
    {
        Vector direction = p2 - p1;

        double angleRadians = angleDegrees * Math.PI / 180;

        Vector leftEdge = RotateVector(direction, -angleRadians);
        Vector rightEdge = RotateVector(direction, angleRadians);
        Point leftPoint = p2 + leftEdge;
        Point rightPoint = p2 + rightEdge;

        Draw(p2, leftPoint);
        Draw(p2, rightPoint);
    }

    private Vector RotateVector(Vector v, double angle)
    {
        double cosTheta = Math.Cos(angle);
        double sinTheta = Math.Sin(angle);
        return new Vector(
            10*(v.X * cosTheta - v.Y * sinTheta),
            10*(v.X * sinTheta + v.Y * cosTheta)
        );
    }

    private void Button_Click(object sender, RoutedEventArgs e)
        {
            switch (((Button)sender).Content)
            {
                case "Clear":
                    clearCanvas();
                    coords.Clear();
                    break;

                case "Toggle Map":
                    image.Visibility = image.Visibility == Visibility.Visible ? Visibility.Hidden : Visibility.Visible;
                    break;
            }
        }

        // Helper Functions
        private void MirrorControlPoints()
        {
            // Mirroring the first and last control point
            var p1 = coords[0];
            var p2 = coords[1];
            coords.Insert(0, new Point(2*p1.X-p2.X, 2* p1.Y - p2.Y));

            p1 = coords[^2];
            p2 = coords[^1];
            
            var last = new Point(2*p2.X - p1.X, 2* p2.Y - p1.Y);
            coords.Add(last);
        }

        private void clearCanvas() => canvas.Children.Clear();

        private void SplineDrawing(Brush colour)
        {
            MirrorControlPoints();
            var splinePoints = cr.CRChain(coords, alpha);
            

            // Remove mirrored control points
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);

            clearCanvas();
            Draw(coords, 2, Brushes.LightGray);
            foreach (var point in coords)
            {
                Draw(point, 6, Brushes.LightGray);
            }
        Draw(splinePoints, 5, colour);
        DrawCone(coords[coords.Count - 2], coords[coords.Count - 1], 45, canvas);

        }

        private void Draw(List<Point> points, int size, Brush colour, bool tracing = true)
        
        {
            if (!tracing)
            {
                foreach (var point in points)
                {
                    Ellipse marker = new Ellipse
                    {
                        Height = size,
                        Width = size,
                        Fill = colour,
                    };
                    Canvas.SetLeft(marker, point.X - marker.Width / 2);
                    Canvas.SetTop(marker, point.Y - marker.Height / 2);
                    canvas.Children.Add(marker);
                }
            }
            else
            {
                for (int i = 0; i < points.Count - 1; i++)
                {
                    Line line = new Line
                    {
                        X1 = points[i].X,
                        Y1 = points[i].Y,
                        X2 = points[i + 1].X,
                        Y2 = points[i + 1].Y,
                        Stroke = colour,
                        StrokeThickness = size
                    };
                    canvas.Children.Add(line);
                }
            }

        }
    
        private void Draw(Point point, int size, Brush colour)
        {
            canvas.Children.Add(
                new Path()
                {
                    Stroke = colour,
                    StrokeThickness = 1,
                    Data = new EllipseGeometry(new Point(point.X, point.Y), size, size)
                }
            );
        }

        private void Draw(Point start, Point end)
        {
            Line line = new Line
            {
                X1 = start.X,
                Y1 = start.Y,
                X2 = end.X,
                Y2 = end.Y,
                Stroke = Brushes.LightGray,
                StrokeThickness = 2
            };
            canvas.Children.Add(line);
        }

        private void slValue_ValueChanged(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            alpha = (float)e.NewValue;
            if (coords.Count > 1){SplineDrawing(Brushes.Orange);}
        }
    }
}
