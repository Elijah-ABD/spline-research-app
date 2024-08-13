using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;

namespace WpfApp1
{
    public partial class MainWindow : System.Windows.Window
    {
        private List<Point> coords = new List<Point>();
        private CatmullRom cr = new CatmullRom();
        private float alpha = 0.5f;
        
        public MainWindow() => InitializeComponent();

        // Event Handlers

        private void MouseClick(object sender, MouseEventArgs e)
        {
            var point = e.GetPosition(this);
            if (coords.Count > 1 && !cr.IsPointInCone(coords[^2], 
                coords[^1], 45, point)) return;
            coords.Add(point);
            Draw(point, 6, Brushes.LightGray);
                  
            if (coords.Count > 1){SplineDrawing(Brushes.Orange);}
        }
        private void sliderChange(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            alpha = (float)e.NewValue;
            if (coords.Count > 1) {SplineDrawing(Brushes.Orange);}
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
        private void Mirror(int index, Point p1, Point p2) => coords.Insert(index, new Point(2 * p1.X - p2.X, 2 * p1.Y - p2.Y));
        private void clearCanvas() => canvas.Children.Clear();

        // Drawing Functions
        private void SplineDrawing(Brush colour)
        {
            Mirror(0, coords[0], coords[1]);
            Mirror(coords.Count, coords[^1], coords[^2]);
            var splinePoints = cr.CRChain(coords, alpha);
            

            // Remove mirrored control points
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);

            clearCanvas();
            Draw(coords, 2, Brushes.LightGray);
            foreach (var point in coords){Draw(point, 6, Brushes.LightGray);}
            
            Draw(splinePoints, 5, colour);
            DrawCone(coords[^2], coords[^1], 45, 2000); // Arbitrary Values

        }

        private void Draw(List<Point> points, int size, Brush colour)
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

        public void DrawCone(Point p1, Point p2, double angleDegrees, double coneLength)
        {
            Vector direction = p2 - p1;
            direction.Normalize();
            double angleRadians = cr.ToRadians(angleDegrees);

            Vector scaledDirection = direction * coneLength;
            Vector leftEdge = cr.RotateVector(scaledDirection, -angleRadians);
            Vector rightEdge = cr.RotateVector(scaledDirection, angleRadians);

            Point leftPoint = p2 + leftEdge;
            Point rightPoint = p2 + rightEdge;

            Draw(p2, leftPoint);
            Draw(p2, rightPoint);
        }
    }
}
