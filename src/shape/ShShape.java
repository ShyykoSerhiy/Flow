package shape;

import draw.Point;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

/**
 * User: shyyko
 * Date: 02.02.14
 * Time: 19:18
 */
public class ShShape extends Shape {
    private static final int AMOUNT_OF_SIDES = 4;
    private int amountOfAdditionalPoints;

    public ShShape(int amountOfAdditionalPoints) {
        assert amountOfAdditionalPoints % AMOUNT_OF_SIDES == 0;
        this.amountOfAdditionalPoints = amountOfAdditionalPoints;

        listOfPoints = new ArrayList<Point>(6);
        listOfPoints.addAll(Arrays.asList(
                new Point(0.8f, 0.2f), new Point(0.8f, 0.8f), new Point(0.2f, 0.8f), new Point(0.2f, 0.2f),
                new Point(0.5f, 0.2f), new Point(0.5f, 0.8f)
        ));
        pointsOfSeparation = new ArrayList<Point>(5);
        pointsOfSeparation.addAll(Arrays.asList(
                listOfPoints.get(0), listOfPoints.get(1), listOfPoints.get(2), listOfPoints.get(3),
                listOfPoints.get(4)
        ));

        kolokPoints = new ArrayList<Point>();

        allPoints = new ArrayList<Point>();
        allPoints.add(listOfPoints.get(0));
        allPoints.addAll(createAdditionalPoints(listOfPoints.get(0), listOfPoints.get(1), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        allPoints.add(listOfPoints.get(1));
        allPoints.addAll(createAdditionalPoints(listOfPoints.get(1), listOfPoints.get(2), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        allPoints.add(listOfPoints.get(2));
        allPoints.addAll(createAdditionalPoints(listOfPoints.get(2), listOfPoints.get(3), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        allPoints.add(listOfPoints.get(3));

        kolokPoints.addAll(createCollocPoints(allPoints));
        //additional side(middle)
        List<Point> middlePoints = new ArrayList<Point>(amountOfAdditionalPoints / AMOUNT_OF_SIDES + 2);
        middlePoints.add(listOfPoints.get(4));
        middlePoints.addAll(createAdditionalPoints(listOfPoints.get(4), listOfPoints.get(5), amountOfAdditionalPoints / AMOUNT_OF_SIDES));
        middlePoints.add(listOfPoints.get(5));
        allPoints.addAll(middlePoints);

        kolokPoints.addAll(createCollocPoints(middlePoints));


        normal = new ArrayList<Point>();
        int kolokPointsPerSide = kolokPoints.size() / AMOUNT_OF_SIDES;
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(1, 0)));
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(0, 1)));
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(-1, 0)));
        normal.addAll(Collections.nCopies(kolokPointsPerSide, new Point(1, 0)));

        assert normal.size() == kolokPoints.size();
    }

    private List<Point> createAdditionalPoints(Point firstPoint, Point secondPoint, int amountOfPoints) {
        List<Point> additionalPoints = new ArrayList<Point>(amountOfPoints);

        Point step = secondPoint.minus(firstPoint).divide(amountOfPoints + 1);
        for (int i = 0; i < amountOfPoints; i++) {
            additionalPoints.add(firstPoint.add(step.multiply(i + 1)));
        }
        return additionalPoints;
    }

    private List<Point> createCollocPoints(List<Point> points) {
        List<Point> kolokPoints = new ArrayList<Point>();
        for (int i = 1; i < points.size(); i++) {
            Point step = points.get(i).minus(points.get(i - 1)).divide(2);
            kolokPoints.add(points.get(i - 1).add(step));
        }
        return kolokPoints;
    }
}
