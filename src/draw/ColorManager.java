package draw;

/**
 * User: shyyko
 * Date: 08.02.14
 * Time: 19:54
 */
public class ColorManager {
    public static final int COLOR_STEPS = 100;

    private double minimalValue;
    private double maximumValue;
    private Color minimumColor;
    private Color maximumColor;

    private double valueStep;
    private Color colorStep;

    public ColorManager(double minimalValue, double maximumValue, Color minimumColor, Color maximumColor) {
        this.minimalValue = minimalValue;
        this.maximumValue = maximumValue;
        this.minimumColor = minimumColor;
        this.maximumColor = maximumColor;

        recalculate();
    }

    public Color getColor(double value) {
        double normalizedValue = value - minimalValue;
        int stepNumberOfValue = (int) Math.ceil(normalizedValue / valueStep);

        return minimumColor.plus(colorStep.multiply(stepNumberOfValue));
    }

    public double getMinimalValue() {
        return minimalValue;
    }

    public void setMinimalValue(double minimalValue) {
        this.minimalValue = minimalValue;
        recalculate();
    }

    public double getMaximumValue() {
        return maximumValue;
    }

    public void setMaximumValue(double maximumValue) {
        this.maximumValue = maximumValue;
        recalculate();
    }

    public Color getMinimumColor() {
        return minimumColor;
    }

    public void setMinimumColor(Color minimumColor) {
        this.minimumColor = minimumColor;
        recalculate();
    }

    public Color getMaximumColor() {
        return maximumColor;
    }

    public void setMaximumColor(Color maximumColor) {
        this.maximumColor = maximumColor;
        recalculate();
    }

    public Color getColorStep() {
        return colorStep;
    }

    public double getValueStep() {
        return valueStep;
    }

    private void recalculate(){
        valueStep = (maximumValue - minimalValue) / COLOR_STEPS;
        colorStep = maximumColor.minus(minimumColor).divide(COLOR_STEPS);
    }
}
