package shape;

import draw.*;
import draw.Color;
import draw.Point;
import processing.core.PApplet;
import solver.Solver;

import java.util.List;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 14:18
 */
public class ShapeDrawer {
    private static final int POINT_SIZE = 10;
    private static final Color RED = new Color(255, 0, 0);
    private static final Color GREEN = new Color(0, 255, 0);
    private static final Color BLUE = new Color(0, 0, 255);

    private Helper helper;
    private Shape shape;
    private PApplet drawer;
    private Solver solver;

    public ShapeDrawer(Helper helper, Shape shape, PApplet drawer, Solver solver) {
        this.helper = helper;
        this.shape = shape;
        this.drawer = drawer;
        this.solver = solver;
    }

    public void draw() {
        drawPoints(shape.getAllPoints(), BLUE);
        drawPoints(shape.getKolokPoints(), RED);
        drawPoints(shape.getListOfPoints(), GREEN);
    }

    private void drawPoints(List<Point> points, Color color) {
        drawer.strokeWeight(POINT_SIZE);
        drawer.stroke(color.getR(), color.getG(), color.getB());
        for (Point point : points) {
            drawPoint(point);
        }
    }

    private void drawPoint(Point coordinatePoint) {
        Point drawPoint = helper.getDrawPoint(coordinatePoint);
        drawer.point((float) drawPoint.getX(), (float) drawPoint.getY());
    }


    private double[][] phiValues;
    private ColorManager colorManager;

    public void computePhi() {
        phiValues = new double[helper.getWidth()][helper.getHeight()];

        double maxPhi = Double.NEGATIVE_INFINITY;
        double minPhi = Double.POSITIVE_INFINITY;
        for (int i = 0; i < phiValues.length; i+= 10) {
            for (int j = 0; j < phiValues[i].length; j+= 10) {
                Point point = helper.getCoordinatePoint(new Point(i, j));
                double phi = solver.getPsi(point);
                if (phi < minPhi) {
                    minPhi = phi;
                }
                if (phi > maxPhi) {
                    maxPhi = phi;
                }
                phiValues[i][j] = phi;
            }
        }
        colorManager = new ColorManager(minPhi, maxPhi, new Color(10, 10, 10), new Color(250, 130, 10));
    }

    public void drawPhi() {
        for (int i = 0; i < phiValues.length; i += 10) {
            for (int j = 0; j < phiValues[i].length; j += 10) {
                Color color = colorManager.getColor(phiValues[i][j]);
                drawer.stroke(color.getR(), color.getG(), color.getB());
                drawer.strokeWeight(POINT_SIZE * 2);
                drawer.point((float) i, (float) j);
            }
        }
    }

    private Point[][] vectorValues;

    public void computeVectors() {
        vectorValues = new Point[helper.getWidth()][helper.getHeight()];
        for (int i = 0; i < vectorValues.length; i += 10) {
            for (int j = 0; j < vectorValues[i].length; j += 10) {
                vectorValues[i][j] = solver.computeV(helper.getCoordinatePoint(new Point(i, j)));
            }
        }
    }

    public void drawVectors() {
        drawer.stroke(10, 10, 10);
        drawer.strokeWeight(2);

        for (int i = 0; i < vectorValues.length; i += 10) {
            for (int j = 0; j < vectorValues[i].length; j += 10) {
                Point vector = vectorValues[i][j].add(new Point(i, j));
                drawer.line((float) i, (float) j, (float) vector.getX(), (float) vector.getY());
            }
        }
    }
}
