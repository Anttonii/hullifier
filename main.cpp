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
    std::pair<int, int> currentPoint;
    std::pair<int, int> nextPoint;
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

double magnitude(std::pair<int, int> a)
{
    return pow(a.first, 2) + pow(a.second, 2);
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
 * Utility function for formatting a point into a string.
 */
std::string formatPoint(std::pair<int, int> p)
{
    return "(" + std::to_string(p.first) + "," + std::to_string(p.second) + ")";
}

/**
 * Point comparison structure that allows the passage of the pivot point into the sorting function.
 */
struct PointCompare
{
    PointCompare(std::pair<int, int> pivot) { this->pivot = pivot; }
    bool operator()(std::pair<int, int> a, std::pair<int, int> b)
    {
        return compareAngles(this->pivot, a, b);
    }

    std::pair<int, int> pivot;
};

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Grahams Scan.
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

        steps.push_back(Step{"Finished.", std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1), points});
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
    std::sort(points.begin() + 1, points.end(), PointCompare(points[0]));

    steps.push_back(Step{"Choosing pivot.", points[0], std::pair(-1, -1), std::vector<std::pair<int, int>>{}});

    // Initially push in the left most and right most angle-wise points w.r.t. point 0.
    std::vector<std::pair<int, int>> output{points[n - 1], points[0]};
    steps.push_back(Step{"Picking rightmost point angle-wise.", points[0], points[n - 1], output});
    output.push_back(points[1]);
    steps.push_back(Step{"Picking leftmost point angle-wise.", points[0], points[1], output});

    int i = 2;
    while (i < n)
    {
        int j = output.size() - 1;
        steps.push_back(Step{"Checking if counter clockwise turn between " + formatPoint(output[j]) + " and " + formatPoint(points[i]) + ".", output[j], points[i], output});

        if (ccw(output[j - 1], output[j], points[i]))
            output.push_back(points[i++]);
        else
            output.pop_back();
    }

    steps.push_back(Step{"Finished.", std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1), output});
    return steps;
}

/**
 * Produces a vector of steps that show the procedure of forming the Convex Hull using Jarvis March.
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

        steps.push_back(Step{"Finished.", std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1), points});
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

    steps.push_back(Step{"Choosing pivot.", points[point0], std::pair(-1, -1), std::vector<std::pair<int, int>>{}});

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

            steps.push_back(Step{"Checking if counter clockwise turn between " + formatPoint(points[currPoint]) + " and " + formatPoint(points[i]) + ".", points[currPoint], points[i], output});
        }

        currPoint = checkPoint;
    }
    // Since current point is now the initial point, add it to the output.
    output.push_back(points[initialPoint]);

    steps.push_back(Step{"Finished.", std::pair<int, int>(-1, -1), std::pair<int, int>(-1, -1), output});
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
void drawLineFromStep(const std::vector<Step> &steps, int index, Color color)
{
    int x1 = steps[index].currentPoint.first;
    int y1 = steps[index].currentPoint.second;
    int x2 = steps[index].nextPoint.first;
    int y2 = steps[index].nextPoint.second;

    if (steps[index].currentPoint != std::pair(-1, -1))
        DrawCircle(x1, y1, 3, color);
    if (steps[index].nextPoint != std::pair(-1, -1))
        DrawCircle(x2, y2, 3, color);
    if (steps[index].currentPoint != std::pair(-1, -1) && steps[index].nextPoint != std::pair(-1, -1))
        DrawLine(x1, y1, x2, y2, color);
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
    if (argc != 2)
    {
        std::cout << "Hullifier only acceps a single command line argument being the path of the input." << std::endl;
        return 0;
    }
    std::string path(argv[1]);

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
        if (IsKeyUp(KEY_SPACE) && keysPressed.count(KEY_SPACE))
        {
            keysPressed.erase(KEY_SPACE);
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
            DrawText("Paused", 10, 500, 18, BLACK);

        drawLineFromStep(steps, currentStep, RED);
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

        DrawText("STEP DURATION:", 510, 155, 20, BLACK);
        GuiSliderBar(Rectangle{555, 180, 105, 20}, "Seconds:", TextFormat("%.2f", stepDuration), &stepDuration, 0.0f, 3.0f);

        DrawText("ALGORITHM:", 535, 95, 20, BLACK);
        if (GuiDropdownBox(Rectangle{505, 120, 180, 20}, "GRAHAM'S SCAN;JARVIS MATCH", &algorithm, algorithmEditMode))
            algorithmEditMode = !algorithmEditMode;

        EndDrawing();
    }

    CloseWindow();

    return 0;
}