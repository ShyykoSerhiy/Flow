package shape;

import draw.Point;

import java.util.List;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 14:18
 */
public abstract class Shape {
    protected List<Point> listOfPoints;
    protected List<Point> allPoints;
    protected List<Point> kolokPoints;
    protected List<Point> normal;
    protected List<Point> pointsOfSeparation;

    public List<Point> getListOfPoints() {
        return listOfPoints;
    }

    public List<Point> getAllPoints() {
        return allPoints;
    }

    public List<Point> getKolokPoints() {
        return kolokPoints;
    }

    public List<Point> getNormal() {
        return normal;
    }

    public List<Point> getPointsOfSeparation() {
        return pointsOfSeparation;
    }
}
