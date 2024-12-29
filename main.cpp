#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <cstdint>
#include <set>
#include <variant>

#include "raylib.h"
#include "raymath.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

/**
 * Represents a point in two-dimensional space.
 */
struct Point2D
{
    int x;
    int y;

    Point2D() : x(0), y(0)
    {
    }

    Point2D(int _x, int _y) : x(_x), y(_y)
    {
    }

    bool operator==(const Point2D &p2)
    {
        return this->x == p2.x && this->y == p2.y;
    }

    bool operator!=(const Point2D &p2)
    {
        return this->x != p2.x || this->y != p2.y;
    }
};

/**
 * Represents a point in three-dimensional space.
 */
struct Point3D
{
    float x;
    float y;
    float z;

    Point3D() : x(0.0f), y(0.0f), z(0.0f)
    {
    }

    Point3D(float _x, float _y, float _z) : x(_x), y(_y), z(_z)
    {
    }

    bool operator==(const Point3D &p2)
    {
        return this->x == p2.x && this->y == p2.y && this->z == p2.z;
    }

    bool operator!=(const Point3D &p2)
    {
        return this->x != p2.x || this->y != p2.y || this->z != p2.z;
    }

    Vector3 toVector() const
    {
        return Vector3{this->x, this->y, this->z};
    }
};

/**
 * Represents a single step in the algorithm when operating in 2-dimensional space.
 */
struct Step2D
{
    std::string stepText;
    std::vector<Point2D> participatingPoints;
    std::vector<Point2D> currentlyAccepted;
};

/**
 * Represents a single step in the algorithm when operating in 2-dimensional space.
 */
struct Step3D
{
    std::string stepText;
    std::vector<Point3D> participatingPoints;
    std::vector<Point3D> currentlyAccepted;
};

/**
 * Utility function for splitting a string.
 */
std::vector<std::string> splitString(std::string &str, char del)
{
    std::stringstream ss;
    std::vector<std::string> output;
    for (auto c : str)
    {
        if (c == del)
        {
            output.push_back(ss.str());
            ss.str("");
        }
        else
        {
            ss << c;
        }
    }
    output.push_back(ss.str());
    return output;
}

/* Variant for the return type of reading a file, since a file can contain either 2D points or 3D points. */
using PointData = std::variant<std::vector<Point2D>, std::vector<Point3D>>;

/**
 * Reads a file from given path and returns all points.
 *
 * A point is a line that has two values separated by a comma.
 */
PointData readFile(std::string &path)
{
    bool is3D = false;
    std::vector<Point2D> points2D;
    std::vector<Point3D> points3D;

    std::ifstream fstream(path);
    if (!fstream.is_open())
        throw std::runtime_error("Failed to open file.");

    std::string line;
    int lineNumber = 0;
    while (std::getline(fstream, line))
    {
        lineNumber++;

        std::vector<std::string> splitted = splitString(line, ',');
        if (splitted.size() < 2 || splitted.size() > 3)
        {
            std::cerr << "Failed to read entry on line: " << lineNumber << std::endl;
            continue;
        }

        // Catch errors when stoi fails.
        try
        {
            if (splitted.size() == 2) // 2D point
            {
                int x = std::stoi(splitted[0]);
                int y = std::stoi(splitted[1]);

                if (x > 500 || x < 0 || y > 500 || y < 0)
                {
                    std::cout << "Invalid point in input, x and y values must be in range [0,500]." << std::endl;
                    continue;
                }
                points2D.push_back(Point2D{x, y});
            }
            else if (splitted.size() == 3) // 3D point
            {
                float x = std::stof(splitted[0]);
                float y = std::stof(splitted[1]);
                float z = std::stof(splitted[2]);

                if (x > 1.0 || x < -1.0 || y > 1.0 || y < -1.0 || z > 1.0 || z < -1.0)
                {
                    std::cout << "Invalid point in input, x, y, z values must be in range [-1,1]." << std::endl;
                    continue;
                }
                points3D.push_back(Point3D{x, y, z});
                is3D = true; // We have detected 3D points
            }
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }

    if (is3D)
    {
        if (points3D.size() == 0)
            throw std::runtime_error("Failed to read input, 0 points loaded, exiting..");

        std::cout << "Loaded " << points3D.size() << " points." << std::endl;
        return points3D;
    }
    else
    {
        if (points2D.size() == 0)
            throw std::runtime_error("Failed to read input, 0 points loaded, exiting..");

        std::cout << "Loaded " << points2D.size() << " points." << std::endl;
        return points2D;
    }
}

/**
 * Get the vector between two points.
 */
Point2D toVector(Point2D a, Point2D b)
{
    return Point2D(b.x - a.x, b.y - a.y);
}

/**
 * Dot product of two vectors
 */
double dot(Point2D a, Point2D b)
{
    return a.x * b.x + a.y * b.y;
}

/**
 * Calculates the magnitude of a given vector.
 */
double magnitude(Point2D a)
{
    return sqrt(pow(a.x, 2) + pow(a.y, 2));
}

/**
 * Cross product of two vectors.
 */
double cross(Point2D a, Point2D b)
{
    return a.x * b.y - a.y * b.x;
}

/**
 * Euclidean distance function.
 */
double distance(Point2D a, Point2D b)
{
    return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

/**
 * Checks whether three points are collinear.
 *
 * Using 10^-12 as a measurement because of floating pointer inaccuracies.
 */
bool pointsCollinear(Point2D p, Point2D a, Point2D b)
{
    return fabs(cross(toVector(p, a), toVector(p, b))) < powf(10, -12);
}

/**
 * Checks whether a turn is counter clockwise or not.
 */
bool ccw(Point2D p, Point2D a, Point2D b)
{
    return cross(toVector(p, a), toVector(p, b)) > 0;
}

/**
 * Function that allows us to compare angles between 3 points.
 *
 * Given a pivot point p and points a and b, find the next counter clockwise turn.
 */
bool compareAngles(Point2D p, Point2D a, Point2D b)
{
    // Special case where points are collinear, choose the closer point.
    if (pointsCollinear(p, a, b))
    {
        return distance(p, a) < distance(p, b);
    }

    Point2D pivotA = toVector(p, a);
    Point2D pivotB = toVector(p, b);
    return (atan2(pivotA.y, pivotA.x) - atan2(pivotB.y, pivotB.x)) < 0;
}

/**
 * Given 3 distinct points, calculate the height of the triangle formed by drawing lines from each to each other.
 *
 * Note that points p and a form the base of the triangle.
 * Note that the height can be negative, this means that the third point has lower Y than the base.
 */
double calculateTriangleHeight(Point2D p, Point2D a, Point2D b)
{
    // Use the shoelace theorem to calculate the area
    // then solve h = 2A/|PQ|
    double area = 0.5 * (p.x * a.y - p.y * a.x + a.x * b.y - a.y * b.x + p.y * b.x - p.x * b.y);
    Point2D vec = toVector(p, a);
    return 2 * area / magnitude(vec);
}

/**
 * Utility function for formatting a point into a string.
 */
std::string formatPoint(Point2D p)
{
    return "(" + std::to_string(p.x) + "," + std::to_string(p.y) + ")";
}

/**
 * Angle comparison structure that allows the passage of the pivot point into the sorting function.
 *
 * Calls the compare angles to sort a set of points in order of their angles to the pivot point.
 */
struct AngleCompare
{
    AngleCompare(Point2D pivot) { this->pivot = pivot; }
    bool operator()(Point2D a, Point2D b)
    {
        return compareAngles(this->pivot, a, b);
    }

    Point2D pivot;
};

/**
 * Side compare determines whether or not a given point is on the right or the left side a line segment.
 */
struct SideCompare
{
    SideCompare(const Point2D &a, const Point2D &b)
    {
        this->pointA = a;
        this->pointB = b;
    }

    bool operator()(const Point2D &p)
    {
        return calculateTriangleHeight(this->pointA, this->pointB, p) > 0;
    }

    Point2D pointA;
    Point2D pointB;
};

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Grahams Scan algorithm.
 *
 * Time complexity for this is n log n.
 */
std::vector<Step2D> getGrahamScanSteps(std::vector<Point2D> &points)
{
    std::vector<Step2D> steps;
    int n = points.size();
    // Edge case check
    if (n <= 3)
    {
        if (!(points[0] == points[n - 1]))
            points.push_back(points[0]);

        steps.push_back(Step2D{"Finished.", std::vector{Point2D(-1, -1), Point2D(-1, -1)}, points});
        return steps;
    }

    // Find the point for which Y is lowest and in case of tie the rightmost w.r.t X.
    // Here lowest technically means highest since Y = 0 is at the top of the window.
    int point0 = 0;
    for (int i = 1; i < n; i++)
    {
        if (points[i].y > points[point0].y || (points[i].y == points[point0].y && points[i].x > points[point0].x))
        {
            point0 = i;
        }
    }

    std::swap(points[0], points[point0]);
    std::sort(points.begin() + 1, points.end(), AngleCompare(points[0]));

    steps.push_back(Step2D{"Choosing pivot.", std::vector{points[0], Point2D(-1, -1)}, std::vector<Point2D>{}});

    // Initially push in the left most and right most angle-wise points w.r.t. point 0.
    std::vector<Point2D> output{points[n - 1], points[0]};
    steps.push_back(Step2D{"Picking rightmost point angle-wise.", std::vector{points[0], points[n - 1]}, output});
    output.push_back(points[1]);
    steps.push_back(Step2D{"Picking leftmost point angle-wise.", std::vector{points[0], points[1]}, output});

    int i = 2;
    while (i < n)
    {
        int j = output.size() - 1;
        steps.push_back(Step2D{"Checking if counter clockwise turn between " + formatPoint(output[j]) + " and " + formatPoint(points[i]) + ".",
                               std::vector{output[j], points[i]}, output});

        if (ccw(output[j - 1], output[j], points[i]))
            output.push_back(points[i++]);
        else
            output.pop_back();
    }

    steps.push_back(Step2D{"Finished.", std::vector{Point2D(-1, -1), Point2D(-1, -1)}, output});
    return steps;
}

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Jarvis March algorithm.
 *
 * Time complexity for this is n*m.
 */
std::vector<Step2D> getJarvisMarchSteps(std::vector<Point2D> &points)
{
    std::vector<Step2D> steps;
    int n = points.size();
    // Edge case check
    if (n <= 3)
    {
        if (!(points[0] == points[n - 1]))
            points.push_back(points[0]);

        steps.push_back(Step2D{"Finished.", std::vector{Point2D(-1, -1), Point2D(-1, -1)}, points});
        return steps;
    }

    // Find the leftmost point to be used as the pivot.
    int point0 = 0;
    for (int i = 1; i < n; i++)
    {
        if (points[i].x < points[point0].x)
        {
            point0 = i;
        }
    }

    steps.push_back(Step2D{"Choosing pivot.", std::vector{points[point0], Point2D(-1, -1)}, std::vector<Point2D>{}});

    std::vector<Point2D> output{};
    int initialPoint = point0;
    output.push_back(points[initialPoint]);

    int currPoint = -1;
    while (currPoint != initialPoint)
    {
        if (currPoint != -1)
            output.push_back(points[currPoint]);

        if (currPoint == -1)
            currPoint = initialPoint;

        int checkPoint = (currPoint + 1) % n;
        for (int i = 0; i < n; i++)
        {
            if (i == checkPoint || i == currPoint)
                continue;

            if (ccw(points[currPoint], points[i], points[checkPoint]))
                checkPoint = i;

            steps.push_back(Step2D{
                "Checking if counter clockwise turn between " + formatPoint(points[currPoint]) + " and " + formatPoint(points[i]) + ".",
                std::vector{points[currPoint], points[i]}, output});
        }

        currPoint = checkPoint;
    }
    // Since current point is now the initial point, add it to the output.
    output.push_back(points[initialPoint]);

    steps.push_back(Step2D{"Finished.", std::vector{Point2D(-1, -1), Point2D(-1, -1)}, output});
    return steps;
}

/**
 * Recursive subroutine of the quickhull algorithm.
 */
void recQuickHull(Point2D p, Point2D q, std::vector<Point2D> points, std::vector<Step2D> &steps, std::vector<Point2D> &output)
{
    int n = points.size();
    if (n == 0)
        return;

    bool right = false;
    double farthestPoint = 0.0;
    int farthestPointIndex = -1;
    for (int i = 0; i < n; i++)
    {
        double distance = calculateTriangleHeight(p, q, points[i]);
        if (fabs(distance) > farthestPoint)
        {
            farthestPoint = fabs(distance);
            farthestPointIndex = i;

            if (distance < 0)
                right = true;
        }
    }

    // Swap places to allow for partitioning.
    std::swap(points[0], points[farthestPointIndex]);
    steps.push_back(Step2D{"Find the most distant point from line segment.", std::vector{p, q, points[0], p}, output});

    // Insert into output.
    if (!right)
    {
        auto outputIndex = std::find(output.begin(), output.end(), p);
        if (outputIndex != output.end())
        {
            outputIndex++;
            output.insert(outputIndex, points[0]);
        }
    }
    else
    {
        auto outputIndex = std::find(output.begin(), output.end(), q);
        if (outputIndex != output.end())
        {
            outputIndex++;
            output.insert(outputIndex, points[0]);
        }
    }

    steps.push_back(Step2D{"Accept the point and continue recursively from the new line segments.", std::vector<Point2D>{}, output});
    auto p1 = std::partition(points.begin() + 1, points.end(), SideCompare(p, points[0]));
    std::vector<Point2D> v1points;
    if (!right)
        v1points = std::vector(points.begin() + 1, p1);
    else
        v1points = std::vector(p1, points.end());

    auto p2 = std::partition(points.begin() + 1, points.end(), SideCompare(q, points[0]));
    std::vector<Point2D> v2points;
    if (!right)
        v2points = std::vector(p2, points.end());
    else
        v2points = std::vector(points.begin() + 1, p2);

    recQuickHull(p, points[0], v1points, steps, output);
    if (right)
        recQuickHull(q, points[0], v2points, steps, output);
    else
        recQuickHull(points[0], q, v2points, steps, output);
}

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Quickhull algorithm.
 *
 * Time complexity for this is n log n.
 */
std::vector<Step2D> getQuickHullSteps(std::vector<Point2D> points)
{
    std::vector<Step2D> steps;
    int n = points.size();

    // Find the leftmost and rightmost points to be used as the pivots.
    int point0 = 0;
    int pointn = 0;

    for (int i = 1; i < n; i++)
    {
        // Leftmost point
        if (points[i].x <= points[point0].x)
        {
            // In case of a tie, make sure the point that is Y-wise lower gets chosen.
            // Since the top edge has Y value of 0, we want the point with higher Y value.
            if (points[i].x == points[point0].x && points[i].y < points[point0].y)
                continue;

            point0 = i;
        }

        // Rightmost point
        if (points[i].x >= points[pointn].x)
        {
            // Same as before but here we choose the point that has the lower Y and is thus the point that is higher.
            if (points[i].x == points[point0].x && points[i].y > points[point0].y)
                continue;

            pointn = i;
        }
    }

    // Swap the point 0 and point n to be the first two points.
    std::swap(points[0], points[point0]);
    std::swap(points[1], points[pointn]);

    // Partition the set of points into two subsets based on which side of the line segment the points are.
    auto partition = std::partition(points.begin() + 2, points.end(), SideCompare(points[0], points[1]));

    std::vector<Point2D> output{points[0], points[1], points[0]};
    steps.push_back(Step2D{"Choosing pivots as leftmost and rightmost points.", std::vector{points[0], points[1]}, output});

    recQuickHull(points[0], points[1], std::vector(points.begin() + 2, partition), steps, output);
    recQuickHull(points[0], points[1], std::vector(partition, points.end()), steps, output);

    steps.push_back(Step2D{"Finished.", std::vector{Point2D(-1, -1), Point2D(-1, -1)}, output});

    return steps;
}

/**
 * Utility function that draws points with radius 3 on the x,y coordinates provided.
 */
void drawDots(std::vector<Point2D> &points, Color color)
{
    for (auto const &p : points)
    {
        int x = p.x;
        int y = p.y;
        DrawCircle(x, y, 3, color);
    }
}

/**
 * Utility function for drawing lines from given steps. Takes into account cases, such as where point is -1, -1 or unavailable.
 */
void drawLinesFromStep(const std::vector<Step2D> &steps, int index, Color color)
{
    int length = steps[index].participatingPoints.size();

    for (int i = 1; i < length; i++)
    {
        auto currentPoint = steps[index].participatingPoints[i - 1];
        auto nextPoint = steps[index].participatingPoints[i];

        int x1 = currentPoint.x;
        int y1 = currentPoint.y;
        int x2 = nextPoint.x;
        int y2 = nextPoint.y;

        if (currentPoint != Point2D(-1, -1))
            DrawCircle(x1, y1, 3, color);
        if (nextPoint != Point2D(-1, -1))
            DrawCircle(x2, y2, 3, color);
        if (currentPoint != Point2D(-1, -1) && nextPoint != Point2D(-1, -1))
            DrawLine(x1, y1, x2, y2, color);
    }
}

/**
 * Draws line between given vector of points with given color.
 *
 * Optionally also draws the points that connect the lines.
 */
void drawLines(const std::vector<Point2D> &points, Color color, bool drawPoints)
{
    for (int i = 1; i < points.size(); i++)
    {
        Point2D p = points[i - 1];
        Point2D p2 = points[i];

        int x1 = p.x;
        int y1 = p.y;
        int x2 = p2.x;
        int y2 = p2.y;
        if (drawPoints)
        {
            DrawCircle(x1, y1, 3, color);
            DrawCircle(x2, y2, 3, color);
        }
        DrawLine(x1, y1, x2, y2, color);
    }
}

/**
 * The main rendering function when operating in 2-dimensions.
 */
void render2D(std::vector<Point2D> &points)
{
    bool running = false;
    auto steps = getGrahamScanSteps(points);
    int currentStep = 0;
    float lastStep = 0;
    float stepDuration = 0.75f;
    std::vector<std::pair<Point2D, Point2D>> lines;

    Rectangle guiBorder = Rectangle{500, 5, 195, 490};

    // Algorithm dropdown selection
    int algorithm = 0;
    int currentAlgorithm = 0;
    bool algorithmEditMode = false;

    // Menu buttons
    bool startToggle = false;
    bool stopToggle = false;
    bool restartToggle = false;

    Rectangle startRectangle = Rectangle{505, 15, 90, 30};
    Rectangle stopRectangle = Rectangle{600, 15, 90, 30};
    Rectangle resetRectangle = Rectangle{505, 50, 185, 30};

    std::set<int> keysPressed;

    while (!WindowShouldClose())
    {
        if (currentAlgorithm != algorithm)
        {
            currentAlgorithm = algorithm;
            running = false;
            lastStep = 0.0f;
            currentStep = 0;
            switch (currentAlgorithm)
            {
            case 0:
                steps = getGrahamScanSteps(points);
                break;
            case 1:
                steps = getJarvisMarchSteps(points);
                break;
            case 2:
                steps = getQuickHullSteps(points);
                break;
            }
        }

        float dt = GetFrameTime();

        if (CheckCollisionPointRec(GetMousePosition(), startRectangle))
        {
            startToggle = true;
            if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
            {
                running = true;
                lastStep = 0.0f;
            }
        }
        else
            startToggle = false;

        if (CheckCollisionPointRec(GetMousePosition(), stopRectangle))
        {
            stopToggle = true;
            if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
            {
                running = false;
                lastStep = 0.0f;
            }
        }
        else
            stopToggle = false;

        if (CheckCollisionPointRec(GetMousePosition(), resetRectangle))
        {
            restartToggle = true;
            if (IsMouseButtonReleased(MOUSE_BUTTON_LEFT))
            {
                running = false;
                currentStep = 0;
                lastStep = 0.0f;
            }
        }
        else
            restartToggle = false;

        if (IsKeyDown(KEY_SPACE) && !keysPressed.count(KEY_SPACE))
        {
            running = !running;
            lastStep = 0.0f;
            keysPressed.insert(KEY_SPACE);
        }

        if (IsKeyDown(KEY_LEFT) && !keysPressed.count(KEY_LEFT))
        {
            lastStep = 0.0f;
            currentStep = std::max(0, currentStep - 1);
            keysPressed.insert(KEY_LEFT);
        }

        if (IsKeyDown(KEY_RIGHT) && !keysPressed.count(KEY_RIGHT))
        {
            lastStep = 0.0f;
            currentStep = std::min((int)steps.size() - 1, currentStep + 1);
            keysPressed.insert(KEY_RIGHT);
        }

        if (IsKeyDown(KEY_F) && !keysPressed.count(KEY_F))
        {
            lastStep = 0.0f;
            currentStep = (int)steps.size() - 1;
            running = false;
            keysPressed.insert(KEY_F);
        }

        if (IsKeyDown(KEY_R) && !keysPressed.count(KEY_R))
        {
            lastStep = 0.0f;
            currentStep = 0;
            running = false;
            keysPressed.insert(KEY_R);
        }

        for (auto const &key : keysPressed)
        {
            if (IsKeyUp(key) && keysPressed.count(key))
                keysPressed.erase(key);
        }

        BeginDrawing();
        ClearBackground(RAYWHITE);

        drawDots(points, BLACK);
        if (running)
        {
            DrawText(std::string("Step: " + steps[currentStep].stepText).c_str(), 10, 500, 18, BLACK);
            if (lastStep >= stepDuration && currentStep < steps.size() - 1)
            {
                currentStep++;
                lastStep = 0;
            }
            lastStep += dt;
        }
        else
        {
            std::string text = currentStep == (int)steps.size() - 1 ? "Finished." : "Paused.";
            DrawText(text.c_str(), 10, 500, 18, BLACK);
        }

        drawLinesFromStep(steps, currentStep, RED);
        drawLines(steps[currentStep].currentlyAccepted, RED, true);

        // Draw GUI
        DrawRectangleRec(startRectangle, startToggle ? LIGHTGRAY : WHITE);
        DrawRectangleLines(startRectangle.x, startRectangle.y, startRectangle.width, startRectangle.height, BLACK);
        DrawRectangleRoundedLinesEx(guiBorder, 0.1f, 0, 1.0f, BLACK);
        DrawText("START", startRectangle.x + 16, startRectangle.y + 8, 18, BLACK);

        DrawRectangleRec(stopRectangle, stopToggle ? LIGHTGRAY : WHITE);
        DrawRectangleLines(stopRectangle.x, stopRectangle.y, stopRectangle.width, stopRectangle.height, BLACK);
        DrawRectangleRoundedLinesEx(guiBorder, 0.1f, 0, 1.0f, BLACK);
        DrawText("STOP", stopRectangle.x + 20, stopRectangle.y + 8, 18, BLACK);

        DrawRectangleRec(resetRectangle, restartToggle ? LIGHTGRAY : WHITE);
        DrawRectangleLines(resetRectangle.x, resetRectangle.y, resetRectangle.width, resetRectangle.height, BLACK);
        DrawRectangleRoundedLinesEx(guiBorder, 0.1f, 0, 1.0f, BLACK);
        DrawText("RESET", resetRectangle.x + 60, resetRectangle.y + 8, 18, BLACK);

        DrawText("STEP DURATION:", 505, 155, 20, BLACK);
        if (!algorithmEditMode)
            GuiSliderBar(Rectangle{555, 180, 105, 20}, "Seconds:", TextFormat("%.2f", stepDuration), &stepDuration, 0.0f, 3.0f);

        DrawText("ALGORITHM:", 535, 95, 20, BLACK);
        if (GuiDropdownBox(Rectangle{505, 120, 180, 20}, "GRAHAM'S SCAN;JARVIS MATCH;QUICK HULL", &algorithm, algorithmEditMode))
            algorithmEditMode = !algorithmEditMode;

        EndDrawing();
    }
}

/**
 * Applies rotation around x and y axis to a vector.
 */
Vector3 applyRotation(Vector3 vec, float rotx, float roty)
{
    Matrix ry = {cos(rotx), 0.0f, sin(rotx), 0.0f,
                 0.0f, 1.0f, 0.0f, 0.0f,
                 -sin(rotx), 0.0f, cos(rotx), 0.0f,
                 0.0f, 0.0f, 0.0f, 1.0f};

    Matrix rx = {1.0f, 0.0f, 0.0f, 0.0f,
                 0.0f, cos(roty), -sin(roty), 0.0f,
                 0.0f, sin(roty), cos(roty), 0.0f,
                 0.0f, 0.0f, 0.0f, 1.0f};

    // Matrix rz = {cos(rotz), -sin(rotz), 0.0f, 0.0f,
    //              sin(rotz), cos(rotz), 0.0f, 0.0f,
    //              0.0f, 0.0f, 1.0f, 0.0f,
    //              0.0f, 0.0f, 0.0f, 1.0f};

    // Matrix rxyz = MatrixMultiply(MatrixMultiply(rz, ry), rx);

    Matrix rxyz = MatrixMultiply(ry, rx);
    return Vector3Transform(vec, rxyz);
}

/**
 * Apply zoom to the the camera when in 3-dimensional space.
 */
Vector3 applyZoom(Vector3 vec, float zoom)
{
    Vector3 min{-1.5f, -1.5f, -3.0f};
    Vector3 max{1.5f, 1.5f, 3.0f};

    return Vector3Clamp(Vector3Scale(vec, zoom), min, max);
}

/**
 * The main rendering function when operating in 3-dimensions.
 */
void render3D(std::vector<Point3D> &points)
{
    Camera3D camera{};
    camera.position = Vector3{0.0f, 1.0f, -1.5f};
    camera.target = Vector3{0.0f, 0.0f, 0.0f};
    camera.up = Vector3{0.0f, 1.0f, 0.0f};
    camera.fovy = 60.0f;
    camera.projection = CAMERA_PERSPECTIVE;

    std::set<int> keysPressed;
    Vector2 lastMouse = Vector2{0.0f, 0.0f};

    while (!WindowShouldClose())
    {
        float dt = GetFrameTime();

        Vector2 mousePos = GetMousePosition();
        Vector2 mouseDelta = Vector2Subtract(mousePos, lastMouse);

        float mouseDeltaX = mouseDelta.x * 0.007f;
        float mouseDeltaY = mouseDelta.y * 0.007f;
        float mouseWheel = -GetMouseWheelMove();

        BeginDrawing();
        ClearBackground(WHITE);
        BeginMode3D(camera);

        for (auto const point : points)
        {
            DrawSphere(point.toVector(), 0.05f, BLACK);
        }

        DrawLine3D(Vector3{-5.0f, 0.0f, 0.0f}, Vector3{5.0f, 0.0f, 0.0f}, BLACK);
        DrawLine3D(Vector3{0.0f, -5.0f, 0.0f}, Vector3{0.0f, 5.0f, 0.0f}, BLACK);
        DrawLine3D(Vector3{0.0f, 0.0f, -5.0f}, Vector3{0.0f, 0.0f, 5.0f}, BLACK);

        EndMode3D();
        EndDrawing();

        if (IsMouseButtonDown(MOUSE_LEFT_BUTTON))
        {
            camera.position = applyRotation(camera.position, mouseDeltaX, mouseDeltaY);
            std::cout << camera.position.x << " " << camera.position.y << " " << camera.position.z << std::endl;
        }

        if (mouseWheel != 0.0f)
        {
            camera.position = applyZoom(camera.position, pow(1.04, mouseWheel));
            std::cout << camera.position.x << " " << camera.position.y << " " << camera.position.z << std::endl;
        }

        for (auto const &key : keysPressed)
        {
            if (IsKeyUp(key) && keysPressed.count(key))
                keysPressed.erase(key);
        }

        lastMouse = mousePos;
    }
}

int main(int argc, char *argv[])
{
    std::string path;
    if (argc == 1)
    {
        std::cout << "Give a path to an input file to start using hullifier." << std::endl;
        std::cin >> path;
        std::cout << "Loading data from path: " << path << std::endl;
    }
    else
        path = std::string(argv[1]);

    if (argc > 2)
    {
        std::cout << "Hullifier only acceps a single command line argument being the path of the input." << std::endl;
        return 0;
    }

    auto points = readFile(path);

    int width = 700;
    int height = 520;
    InitWindow(width, height, "Hullifier");

    if (std::holds_alternative<std::vector<Point2D>>(points))
    {
        auto points2D = std::get<std::vector<Point2D>>(points);
        render2D(points2D);
    }
    else if (std::holds_alternative<std::vector<Point3D>>(points))
    {
        auto points3D = std::get<std::vector<Point3D>>(points);
        render3D(points3D);
    }

    CloseWindow();

    return 0;
}