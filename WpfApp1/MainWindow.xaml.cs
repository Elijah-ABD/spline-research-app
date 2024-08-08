﻿using System.Windows;
using System.Windows.Controls;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Shapes;
using MathNet.Numerics;

namespace WpfApp1
{
    public partial class MainWindow : System.Windows.Window
    {
        private List<Point> coords = new List<Point>();
        private CatmullRom cr = new CatmullRom();
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
                MirrorControlPoints();
                clearCanvas();
                SplineDrawing(Brushes.Yellow, 0.5f);
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
           clearCanvas();
           coords.Clear();
        }

        // Helper Functions
        private void MirrorControlPoints()
        {
            var p1 = coords[0];
            var p2 = coords[1];
            coords.Insert(0, new Point(2*p1.X-p2.X, 2* p1.Y - p2.Y));

            p1 = coords[^2];
            p2 = coords[^1];
            
            var last = new Point(2*p2.X - p1.X, 2* p2.Y - p1.Y);
            coords.Add(last);
        }

        private void clearCanvas()
        {
                canvas.Children.Clear();
        }

        private void SplineDrawing(Brush colour, float alpha = 0)
        {
            var splinePoints = cr.CRChain(coords, alpha);
            coords.RemoveAt(0);
            coords.RemoveAt(coords.Count - 1);
            Draw(coords, 2, Brushes.LightGray);
            foreach (var point in coords)
            {
                Draw(point, 2, colour);
            }
            Draw(splinePoints, 2, colour);

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
                    Stroke = Brushes.LightGray,
                    StrokeThickness = 1,
                    Data = new EllipseGeometry(new Point(point.X, point.Y), 6, 6)
                }
                );
        } 
    }
}
