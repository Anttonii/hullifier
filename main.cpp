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

#include "raylib.h"

#define RAYGUI_IMPLEMENTATION
#include "raygui.h"

/**
 * Represents a single step in the algorithm.
 */
struct Step
{
    std::string stepText;
    std::vector<std::pair<int, int>> participatingPoints;
    std::vector<std::pair<int, int>> currentlyAccepted;
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

/**
 * Reads a file from given path and returns all points.
 *
 * A point is a line that has two values separated by a comma.
 */
std::vector<std::pair<int, int>> readFile(std::string &path)
{
    std::vector<std::pair<int, int>> points;
    std::ifstream fstream(path);
    if (!fstream.is_open())
    {
        std::cerr << "Failed to open file." << std::endl;
        return points;
    }

    std::string line;
    int lineNumber = 0;
    while (std::getline(fstream, line))
    {
        lineNumber++;

        std::vector<std::string> splitted = splitString(line, ',');
        if (splitted.size() != 2)
        {
            std::cerr << "Failed to read entry on line: " << lineNumber << std::endl;
            continue;
        }

        // Catch errors when stoi fails.
        try
        {
            int x = std::stoi(splitted[0]);
            int y = std::stoi(splitted[1]);

            if (x > 500 || x < 0 || y > 500 || y < 0)
                std::cout << "Invalid point in input, x and y values must be in ranged [0,500]." << std::endl;

            std::pair<int, int> point{x, y};
            points.push_back(point);
        }
        catch (const std::exception &e)
        {
            std::cerr << e.what() << std::endl;
        }
    }

    std::cout << "Loaded " << points.size() << " points." << std::endl;

    return points;
}

/**
 * Get the vector between two points.
 */
std::pair<int, int> toVector(std::pair<int, int> a, std::pair<int, int> b)
{
    return std::pair(b.first - a.first, b.second - a.second);
}

/**
 * Dot product of two vectors
 */
double dot(std::pair<int, int> a, std::pair<int, int> b)
{
    return a.first * b.first + a.second * b.second;
}

/**
 * Calculates the magnitude of a given vector.
 */
double magnitude(std::pair<int, int> a)
{
    return sqrt(pow(a.first, 2) + pow(a.second, 2));
}

/**
 * Cross product of two vectors.
 */
double cross(std::pair<int, int> a, std::pair<int, int> b)
{
    return a.first * b.second - a.second * b.first;
}

/**
 * Euclidean distance function.
 */
double distance(std::pair<int, int> a, std::pair<int, int> b)
{
    return sqrt(pow(a.first - b.first, 2) + pow(a.second - b.second, 2));
}

/**
 * Checks whether three points are collinear.
 *
 * Using 10^-12 as a measurement because of floating pointer inaccuracies.
 */
bool pointsCollinear(std::pair<int, int> p, std::pair<int, int> a, std::pair<int, int> b)
{
    return fabs(cross(toVector(p, a), toVector(p, b))) < powf(10, -12);
}

/**
 * Checks whether a turn is counter clockwise or not.
 */
bool ccw(std::pair<int, int> p, std::pair<int, int> a, std::pair<int, int> b)
{
    return cross(toVector(p, a), toVector(p, b)) > 0;
}

/**
 * Function that allows us to compare angles between 3 points.
 *
 * Given a pivot point p and points a and b, find the next counter clockwise turn.
 */
bool compareAngles(std::pair<int, int> p, std::pair<int, int> a, std::pair<int, int> b)
{
    // Special case where points are collinear, choose the closer point.
    if (pointsCollinear(p, a, b))
    {
        return distance(p, a) < distance(p, b);
    }

    std::pair<int, int> pivotA = toVector(p, a);
    std::pair<int, int> pivotB = toVector(p, b);
    return (atan2(pivotA.second, pivotA.first) - atan2(pivotB.second, pivotB.first)) < 0;
}

/**
 * Given 3 distinct points, calculate the height of the triangle formed by drawing lines from each to each other.
 *
 * Note that points p and a form the base of the triangle.
 * Note that the height can be negative, this means that the third point has lower Y than the base.
 */
double calculateTriangleHeight(std::pair<int, int> p, std::pair<int, int> a, std::pair<int, int> b)
{
    // Use the shoelace theorem to calculate the area
    // then solve h = 2A/|PQ|
    double area = 0.5 * (p.first * a.second - p.second * a.first + a.first * b.second - a.second * b.first + p.second * b.first - p.first * b.second);
    std::pair<int, int> vec = toVector(p, a);
    return 2 * area / magnitude(vec);
}

/**
 * Utility function for formatting a point into a string.
 */
std::string formatPoint(std::pair<int, int> p)
{
    return "(" + std::to_string(p.first) + "," + std::to_string(p.second) + ")";
}

/**
 * Angle comparison structure that allows the passage of the pivot point into the sorting function.
 *
 * Calls the compare angles to sort a set of points in order of their angles to the pivot point.
 */
struct AngleCompare
{
    AngleCompare(std::pair<int, int> pivot) { this->pivot = pivot; }
    bool operator()(std::pair<int, int> a, std::pair<int, int> b)
    {
        return compareAngles(this->pivot, a, b);
    }

    std::pair<int, int> pivot;
};

/**
 * Side compare determines whether or not a given point is on the right or the left side a line segment.
 */
struct SideCompare
{
    SideCompare(const std::pair<int, int> &a, const std::pair<int, int> &b)
    {
        this->pointA = a;
        this->pointB = b;
    }

    bool operator()(const std::pair<int, int> &p)
    {
        return calculateTriangleHeight(this->pointA, this->pointB, p) > 0;
    }

    std::pair<int, int> pointA;
    std::pair<int, int> pointB;
};

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Grahams Scan algorithm.
 *
 * Time complexity for this is n log n.
 */
std::vector<Step> getGrahamScanSteps(std::vector<std::pair<int, int>> &points)
{
    std::vector<Step> steps;
    int n = points.size();
    // Edge case check
    if (n <= 3)
    {
        if (!(points[0] == points[n - 1]))
            points.push_back(points[0]);

        steps.push_back(Step{"Finished.", std::vector{std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1)}, points});
        return steps;
    }

    // Find the point for which Y is lowest and in case of tie the rightmost w.r.t X.
    // Here lowest technically means highest since Y = 0 is at the top of the window.
    int point0 = 0;
    for (int i = 1; i < n; i++)
    {
        if (points[i].second > points[point0].second || (points[i].second == points[point0].second && points[i].first > points[point0].first))
        {
            point0 = i;
        }
    }

    std::swap(points[0], points[point0]);
    std::sort(points.begin() + 1, points.end(), AngleCompare(points[0]));

    steps.push_back(Step{"Choosing pivot.", std::vector{points[0], std::pair(-1, -1)}, std::vector<std::pair<int, int>>{}});

    // Initially push in the left most and right most angle-wise points w.r.t. point 0.
    std::vector<std::pair<int, int>> output{points[n - 1], points[0]};
    steps.push_back(Step{"Picking rightmost point angle-wise.", std::vector{points[0], points[n - 1]}, output});
    output.push_back(points[1]);
    steps.push_back(Step{"Picking leftmost point angle-wise.", std::vector{points[0], points[1]}, output});

    int i = 2;
    while (i < n)
    {
        int j = output.size() - 1;
        steps.push_back(Step{"Checking if counter clockwise turn between " + formatPoint(output[j]) + " and " + formatPoint(points[i]) + ".",
                             std::vector{output[j], points[i]}, output});

        if (ccw(output[j - 1], output[j], points[i]))
            output.push_back(points[i++]);
        else
            output.pop_back();
    }

    steps.push_back(Step{"Finished.", std::vector{std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1)}, output});
    return steps;
}

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Jarvis March algorithm.
 *
 * Time complexity for this is n*m.
 */
std::vector<Step> getJarvisMarchSteps(std::vector<std::pair<int, int>> &points)
{
    std::vector<Step> steps;
    int n = points.size();
    // Edge case check
    if (n <= 3)
    {
        if (!(points[0] == points[n - 1]))
            points.push_back(points[0]);

        steps.push_back(Step{"Finished.", std::vector{std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1)}, points});
        return steps;
    }

    // Find the leftmost point to be used as the pivot.
    int point0 = 0;
    for (int i = 1; i < n; i++)
    {
        if (points[i].first < points[point0].first)
        {
            point0 = i;
        }
    }

    steps.push_back(Step{"Choosing pivot.", std::vector{points[point0], std::pair(-1, -1)}, std::vector<std::pair<int, int>>{}});

    std::vector<std::pair<int, int>> output{};
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

            steps.push_back(Step{
                "Checking if counter clockwise turn between " + formatPoint(points[currPoint]) + " and " + formatPoint(points[i]) + ".",
                std::vector{points[currPoint], points[i]}, output});
        }

        currPoint = checkPoint;
    }
    // Since current point is now the initial point, add it to the output.
    output.push_back(points[initialPoint]);

    steps.push_back(Step{"Finished.", std::vector{std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1)}, output});
    return steps;
}

/**
 * Recursive subroutine of the quickhull algorithm.
 */
void recQuickHull(std::pair<int, int> p, std::pair<int, int> q, std::vector<std::pair<int, int>> points, std::vector<Step> &steps, std::vector<std::pair<int, int>> &output)
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
    steps.push_back(Step{"Find the most distant point from line segment.", std::vector{p, q, points[0], p}, output});

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

    steps.push_back(Step{"Accept the point and continue recursively from the new line segments.", std::vector<std::pair<int, int>>{}, output});
    auto p1 = std::partition(points.begin() + 1, points.end(), SideCompare(p, points[0]));
    std::vector<std::pair<int, int>> v1points;
    if (!right)
        v1points = std::vector(points.begin() + 1, p1);
    else
        v1points = std::vector(p1, points.end());

    auto p2 = std::partition(points.begin() + 1, points.end(), SideCompare(q, points[0]));
    std::vector<std::pair<int, int>> v2points;
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
std::vector<Step> getQuickHullSteps(std::vector<std::pair<int, int>> points)
{
    std::vector<Step> steps;
    int n = points.size();

    // Find the leftmost and rightmost points to be used as the pivots.
    int point0 = 0;
    int pointn = 0;

    for (int i = 1; i < n; i++)
    {
        // Leftmost point
        if (points[i].first <= points[point0].first)
        {
            // In case of a tie, make sure the point that is Y-wise lower gets chosen.
            // Since the top edge has Y value of 0, we want the point with higher Y value.
            if (points[i].first == points[point0].first && points[i].second < points[point0].second)
                continue;

            point0 = i;
        }

        // Rightmost point
        if (points[i].first >= points[pointn].first)
        {
            // Same as before but here we choose the point that has the lower Y and is thus the point that is higher.
            if (points[i].first == points[point0].first && points[i].second > points[point0].second)
                continue;

            pointn = i;
        }
    }

    // Swap the point 0 and point n to be the first two points.
    std::swap(points[0], points[point0]);
    std::swap(points[1], points[pointn]);

    // Partition the set of points into two subsets based on which side of the line segment the points are.
    auto partition = std::partition(points.begin() + 2, points.end(), SideCompare(points[0], points[1]));

    std::vector<std::pair<int, int>> output{points[0], points[1], points[0]};
    steps.push_back(Step{"Choosing pivots as leftmost and rightmost points.", std::vector{points[0], points[1]}, output});

    recQuickHull(points[0], points[1], std::vector(points.begin() + 2, partition), steps, output);
    recQuickHull(points[0], points[1], std::vector(partition, points.end()), steps, output);

    steps.push_back(Step{"Finished.", std::vector{std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1)}, output});

    return steps;
}

/**
 * Utility function that draws points with radius 3 on the x,y coordinates provided.
 */
void drawDots(std::vector<std::pair<int, int>> &points, Color color)
{
    for (auto const &p : points)
    {
        int x = p.first;
        int y = p.second;
        DrawCircle(x, y, 3, color);
    }
}

/**
 * Utility function for drawing lines from given steps. Takes into account cases, such as where point is -1, -1 or unavailable.
 */
void drawLinesFromStep(const std::vector<Step> &steps, int index, Color color)
{
    int length = steps[index].participatingPoints.size();

    for (int i = 1; i < length; i++)
    {
        auto currentPoint = steps[index].participatingPoints[i - 1];
        auto nextPoint = steps[index].participatingPoints[i];

        int x1 = currentPoint.first;
        int y1 = currentPoint.second;
        int x2 = nextPoint.first;
        int y2 = nextPoint.second;

        if (currentPoint != std::pair(-1, -1))
            DrawCircle(x1, y1, 3, color);
        if (nextPoint != std::pair(-1, -1))
            DrawCircle(x2, y2, 3, color);
        if (currentPoint != std::pair(-1, -1) && nextPoint != std::pair(-1, -1))
            DrawLine(x1, y1, x2, y2, color);
    }
}

/**
 * Draws line between given vector of points with given color.
 *
 * Optionally also draws the points that connect the lines.
 */
void drawLines(const std::vector<std::pair<int, int>> &points, Color color, bool drawPoints)
{
    for (int i = 1; i < points.size(); i++)
    {
        std::pair<int, int> p = points[i - 1];
        std::pair<int, int> p2 = points[i];

        int x1 = p.first;
        int y1 = p.second;
        int x2 = p2.first;
        int y2 = p2.second;
        if (drawPoints)
        {
            DrawCircle(x1, y1, 3, color);
            DrawCircle(x2, y2, 3, color);
        }
        DrawLine(x1, y1, x2, y2, color);
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
    if (points.size() == 0)
    {
        std::cout << "Failed to read input, 0 points loaded, exiting.." << std::endl;
        return 0;
    }

    auto steps = getGrahamScanSteps(points);
    int currentStep = 0;
    float lastStep = 0;
    float stepDuration = 0.75f;
    std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> lines;

    int width = 700;
    int height = 520;
    InitWindow(width, height, "Hullifier");
    bool running = false;

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

    CloseWindow();

    return 0;
}