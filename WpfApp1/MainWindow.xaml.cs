using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
namespace WpfApp1
{
    public partial class MainWindow : Window
    {
        private List<Point> coords = new List<Point>();
        private List<Path> paths = new List<Path>();
        private CatmullRom cr = new CatmullRom();
        private Drawer dr;
        private float alpha = 0.5f;
        private bool draw = true;
        private Nullable<Point> dragStart = null;

        public MainWindow() {InitializeComponent(); dr = new Drawer(canvas, coords, paths);}

        // Event Handlers
        private void MouseClick(object sender, MouseEventArgs e)
        {
            if (draw)
            {
                var point = e.GetPosition(this);
                if (coords.Count > 1 && !cr.IsPointInCone(coords[^2],
                coords[^1], 45, point)) return;
                coords.Add(point);
                paths.Add(dr.Draw(point, 6, Brushes.LightGray));
                if (coords.Count > 1) { dr.SplineDrawing(Brushes.Orange, alpha, draw); }
            }
        }
        private void sliderChange(object sender, RoutedPropertyChangedEventArgs<double> e)
        {
            alpha = (float)e.NewValue;
            if (coords.Count > 1) {dr.SplineDrawing(Brushes.Orange, alpha, draw);}
        }


        private void ClearButtonClick(object sender, RoutedEventArgs e) { dr.clearCanvas(); coords.Clear(); paths.Clear(); }

        private void DrawButtonClick(object sender, RoutedEventArgs e) 
        { 
            draw = !draw;
            foreach (var path in paths)
            {
                Drag(path, draw);
            }

            if (!draw && canvas.Children.Count > 0 && canvas.Children[^1] is Line)
            {
                canvas.Children.RemoveAt(canvas.Children.Count - 1);
                canvas.Children.RemoveAt(canvas.Children.Count - 1);
            }

            else if (draw && coords.Count > 1)
            {
               dr.Draw(coords[^2], coords[^1], 45, 2000); 
            }
        }
        
        private void MapButtonClick(object sender, RoutedEventArgs e) => image.Visibility = image.Visibility
                                                                      == Visibility.Visible 
                                                                      ? Visibility.Hidden : Visibility.Visible;
        
        private void UndoButtonClick(object sender, RoutedEventArgs e)
        {
            switch (coords.Count)
            {
                case > 2:
                    coords.RemoveAt(coords.Count - 1);
                    paths.RemoveAt(paths.Count - 1);
                    dr.SplineDrawing(Brushes.Orange, alpha, draw);
                    break;

                case 2:
                    coords.RemoveAt(coords.Count - 1);
                    dr.clearCanvas();
                    dr.Draw(coords[0], 6, Brushes.LightGray);
                    break;

                default:
                    dr.clearCanvas();
                    coords.Clear();
                    paths.Clear();
                    break;
            }
        }

        private void Down(object sender, MouseEventArgs e)
        {
            var element = (UIElement)sender;
            dragStart = e.GetPosition(element);
            element.CaptureMouse();
        }
        
        private void Up(object sender, MouseEventArgs e)
        {
            var element = (UIElement)sender;
            dragStart = null;
            element.ReleaseMouseCapture();
        }

        private void Move(object sender, MouseEventArgs e)
        {
            if (dragStart != null && e.LeftButton == MouseButtonState.Pressed)
            {
                var element = (UIElement)sender;
                var p2 = e.GetPosition(this);
                Canvas.SetLeft(element, p2.X - dragStart.Value.X);
                Canvas.SetTop(element, p2.Y - dragStart.Value.Y);

                coords[paths.IndexOf((Path)element)] = p2;
                if (coords.Count > 1)
                {
                    if (cr.IsSplineValid(coords, 45))
                    {
                        dr.SplineDrawing(Brushes.Orange, alpha, draw);
                    }
                    else
                    {
                        dr.SplineDrawing(Brushes.Red, alpha, draw);
                    }
                }
            }
        }
        
        private void Drag(UIElement element, bool draw)
        {
            if (draw)
            {
                element.MouseLeftButtonDown -= Down; 
                element.MouseLeftButtonUp -= Up; 
                element.MouseMove -= Move;
            }
            else
            {
                element.MouseLeftButtonDown += Down;
                element.MouseLeftButtonUp += Up;
                element.MouseMove += Move;
            }
        }
    }
}
