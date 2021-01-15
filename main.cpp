//Algorithms taken from:
//Grahm Scan: https://www.geeksforgeeks.org/convex-hull-set-2-graham-scan/
//Jarvis March: https://www.geeksforgeeks.org/convex-hull-set-1-jarviss-algorithm-or-wrapping/
//Quickhull: https://www.geeksforgeeks.org/quickhull-algorithm-convex-hull/ (edited to work with other implementations)

#include <string>
#include <stdlib.h>
#include <vector>
#include <stack>
#include <set>
#include <iostream>
#include <fstream>
#include <time.h>
#include <bits/stdc++.h>
using std::cout;
using std::vector;
using std::endl;
using std::pair;

struct Point
{
	int x, y;
};

////////////////////////////////////////////////////////////////////////////////////////
// Jarvis March Setup																  //
////////////////////////////////////////////////////////////////////////////////////////

// To find orientation of ordered triplet (p, q, r).
// The function returns following values
// 0 --> p, q and r are colinear
// 1 --> Clockwise
// 2 --> Counterclockwise
int orientation(Point p, Point q, Point r)
{
	int val = (q.y - p.y) * (r.x - q.x) -
		(q.x - p.x) * (r.y - q.y);

	if (val == 0) return 0; // colinear
	return (val > 0) ? 1 : 2; // clock or counterclock wise
}

// Prints convex hull of a set of n points. 
void jarvisHull(Point points[], int n)
{
	// There must be at least 3 points 
	if (n < 3) return;

	// Initialize Result 
	vector<Point> hull;

	// Find the leftmost point 
	int l = 0;
	for (int i = 1; i < n; i++)
		if (points[i].x < points[l].x)
			l = i;

	// Start from leftmost point, keep moving counterclockwise 
	// until reach the start point again.  This loop runs O(h) 
	// times where h is number of points in result or output. 
	int p = l, q;
	do
	{
		// Add current point to result 
		hull.push_back(points[p]);

		// Search for a point 'q' such that orientation(p, x, 
		// q) is counterclockwise for all points 'x'. The idea 
		// is to keep track of last visited most counterclock- 
		// wise point in q. If any point 'i' is more counterclock- 
		// wise than q, then update q. 
		q = (p + 1) % n;
		for (int i = 0; i < n; i++)
		{
			// If i is more counterclockwise than current q, then 
			// update q 
			if (orientation(points[p], points[i], points[q]) == 2)
				q = i;
		}

		// Now q is the most counterclockwise with respect to p 
		// Set p as q for next iteration, so that q is added to 
		// result 'hull' 
		p = q;

	} while (p != l);  // While we don't come to first point 

	// Output hull to file
	std::ofstream ofile("hull_J.txt");
	for (int i = 0; i < hull.size(); i++) {
		int x, y;
		x = hull[i].x;
		y = hull[i].y;
		ofile << x << "\t" << y << "\n";
	}
}

////////////////////////////////////////////////////////////////////////////////////////
// Grahm Scan Setup																	  //
////////////////////////////////////////////////////////////////////////////////////////

// A global point needed for  sorting points with reference
// to  the first point Used in compare function of qsort()
Point p0;

// A utility function to find next to top in a stack
Point nextToTop(std::stack<Point>& S)
{
	Point p = S.top();
	S.pop();
	Point res = S.top();
	S.push(p);
	return res;
}

// A utility function to swap two points
void swap(Point& p1, Point& p2)
{
	Point temp = p1;
	p1 = p2;
	p2 = temp;
}

// A utility function to return square of distance
// between p1 and p2
int distSq(Point p1, Point p2)
{
	return (p1.x - p2.x) * (p1.x - p2.x) +
		(p1.y - p2.y) * (p1.y - p2.y);
}

// A function used by library function qsort() to sort an array of
// points with respect to the first point
int compare(const void* vp1, const void* vp2)
{
	Point* p1 = (Point*)vp1;
	Point* p2 = (Point*)vp2;

	// Find orientation
	int o = orientation(p0, *p1, *p2);
	if (o == 0)
		return (distSq(p0, *p2) >= distSq(p0, *p1)) ? -1 : 1;

	return (o == 2) ? -1 : 1;
}

// Prints convex hull of a set of n points.
void grahmHull(Point points[], int n)
{
	// Find the bottommost point
	int ymin = points[0].y, min = 0;
	for (int i = 1; i < n; i++)
	{
		int y = points[i].y;

		// Pick the bottom-most or chose the left
		// most point in case of tie
		if ((y < ymin) || (ymin == y &&
			points[i].x < points[min].x))
			ymin = points[i].y, min = i;
	}

	// Place the bottom-most point at first position
	swap(points[0], points[min]);

	// Sort n-1 points with respect to the first point.
	// A point p1 comes before p2 in sorted output if p2
	// has larger polar angle (in counterclockwise
	// direction) than p1
	p0 = points[0];
	qsort(&points[1], n - 1, sizeof(Point), compare);

	// If two or more points make same angle with p0,
	// Remove all but the one that is farthest from p0
	// Remember that, in above sorting, our criteria was
	// to keep the farthest point at the end when more than
	// one points have same angle.
	int m = 1; // Initialize size of modified array
	for (int i = 1; i < n; i++)
	{
		// Keep removing i while angle of i and i+1 is same
		// with respect to p0
		while (i < n - 1 && orientation(p0, points[i],
			points[i + 1]) == 0)
			i++;


		points[m] = points[i];
		m++;  // Update size of modified array
	}

	// If modified array of points has less than 3 points,
	// convex hull is not possible
	if (m < 3) return;

	// Create an empty stack and push first three points
	// to it.
	std::stack<Point> S;
	S.push(points[0]);
	S.push(points[1]);
	S.push(points[2]);

	// Process remaining n-3 points
	for (int i = 3; i < m; i++)
	{
		// Keep removing top while the angle formed by
		// points next-to-top, top, and points[i] makes
		// a non-left turn
		while (orientation(nextToTop(S), S.top(), points[i]) != 2)
			S.pop();
		S.push(points[i]);
	}

	std::ofstream ofile("hull_G.txt");
	// Now stack has the output points, output contents of stack
	while (!S.empty())
	{
		Point p = S.top();
		ofile << p.x << "\t" << p.y << endl;
		S.pop();
	}
}

////////////////////////////////////////////////////////////////////////////////////////
// Quickhull Setup																	  //
////////////////////////////////////////////////////////////////////////////////////////

// iPair is integer pairs 
#define iPair pair<int, int> 

// Stores the result (points of convex hull) 
std::set<iPair> hull;

// Returns the side of point p with respect to line 
// joining points p1 and p2. 
int findSide(iPair p1, iPair p2, iPair p)
{
	int val = (p.second - p1.second) * (p2.first - p1.first) -
		(p2.second - p1.second) * (p.first - p1.first);

	if (val > 0)
		return 1;
	if (val < 0)
		return -1;
	return 0;
}

// returns a value proportional to the distance 
// between the point p and the line joining the 
// points p1 and p2 
int lineDist(iPair p1, iPair p2, iPair p)
{
	return abs((p.second - p1.second) * (p2.first - p1.first) -
		(p2.second - p1.second) * (p.first - p1.first));
}

// End points of line L are p1 and p2.  side can have value 
// 1 or -1 specifying each of the parts made by the line L 
void quickHull(iPair a[], int n, iPair p1, iPair p2, int side)
{
	int ind = -1;
	int max_dist = 0;

	// finding the point with maximum distance 
	// from L and also on the specified side of L. 
	for (int i = 0; i < n; i++)
	{
		int temp = lineDist(p1, p2, a[i]);
		if (findSide(p1, p2, a[i]) == side && temp > max_dist)
		{
			ind = i;
			max_dist = temp;
		}
	}

	// If no point is found, add the end points 
	// of L to the convex hull. 
	if (ind == -1)
	{
		hull.insert(p1);
		hull.insert(p2);
		return;
	}

	// Recur for the two parts divided by a[ind] 
	quickHull(a, n, a[ind], p1, -findSide(a[ind], p1, p2));
	quickHull(a, n, a[ind], p2, -findSide(a[ind], p2, p1));
}

void printHull(iPair a[], int n)
{
	// a[i].second -> y-coordinate of the ith point 
	if (n < 3)
	{
		cout << "Convex hull not possible\n";
		return;
	}

	// Finding the point with minimum and 
	// maximum x-coordinate 
	int min_x = 0, max_x = 0;
	for (int i = 1; i < n; i++)
	{
		if (a[i].first < a[min_x].first)
			min_x = i;
		if (a[i].first > a[max_x].first)
			max_x = i;
	}

	// Recursively find convex hull points on 
	// one side of line joining a[min_x] and 
	// a[max_x] 
	quickHull(a, n, a[min_x], a[max_x], 1);

	// Recursively find convex hull points on 
	// other side of line joining a[min_x] and 
	// a[max_x] 
	quickHull(a, n, a[min_x], a[max_x], -1);

	// Output hull to a file
	std::ofstream ofile("hull_Q.txt");
	while (!hull.empty())
	{
		ofile << (*hull.begin()).first << "\t"
			<< (*hull.begin()).second << "\n";
		hull.erase(hull.begin());
	}
}

int main(int argc, char * argv[]){

	if (argc < 3)
		cout << "wrong format! should be \"a.out algType dataFile";
	else {
		//set the algorithm type to be the argument brought in at argv[1]
		std::string algType = argv[1];
		//set the filename for the data to be the argument brought in at argv[2]
		std::string dataFilename = argv[2];
		//output the filename
		cout << "dataFilename = " << dataFilename << endl;
		//========================= Find File Size ======================================
		//could be done more gracefully, but it works.
		//open the data file
		std::ifstream ifile(dataFilename);
		//temporary integers
		int p1, p2;
		//counter value initialized to 0
		int i = 0;
		//while values are still being brought in, increment i
		while (ifile >> p1 >> p2) {
			i++;
		}
		//i is now the number of points
		int numberPoints = i;
		//output the number of points
		cout << "number of points in the file: " << numberPoints << endl;

		//========================= Bring in Points from File ===========================
		//point array to store all the points
		Point points[numberPoints];
		//reset the counter
		i = 0;
		//open the data file again
		std::ifstream readfile(dataFilename);
		//temporary integers
		int x, y;
		//while values are still being brought in...
		while (readfile >> x >> y) {
			//create a temporary point
			Point temp;
			//set the temporary x value to the x value from the file
			temp.x = x;
			//set the temporary y value to the y value from the file
			temp.y = y;
			//set the point at the current count to be the temporary point
			points[i] = temp;
			//increment the counter
			i++;
		}

		//========================= Run Hull Algorithm ==================================
		//if the algorithm type is Grahm scan...
		if (algType == "G") {
			//let the user know that it's grahm scan (this made it easier to quickly look at all the times)
			cout << "\nGrahm Scan\n";
			//create the clock_t value time
			clock_t time;
			//set the clock
			time = clock();
			//generate the hull using Grahm Scan
			grahmHull(points, numberPoints);
			//set time to be the current time minus the start time
			time = clock() - time;
			//output the time in ticks, then seconds
			cout << "completed in " << time << " ticks (" << ((float)time) / CLOCKS_PER_SEC << ")\n";
		}
		//if the algorithm type is Jarvis March...
		else if (algType == "J") {
			//let the user know that it's Jarvis March (this made it easier to quickly look at all the times)
			cout << "\nJarvis March\n";
			//create the clock_t value time
			clock_t time;
			//set the clock
			time = clock();
			//generate the hull using Jarvis March
			jarvisHull(points, numberPoints);
			//set time to be the current time minus the start time
			time = clock() - time;
			//output the time in ticks, then seconds
			cout << "completed in " << time << " ticks (" << ((float)time) / CLOCKS_PER_SEC << ")\n";
		}
		//if the algorithm type is Quickhull...
		else if (algType == "Q"){
			//let the user know that it's Quickhull (this made it easier to quickly look at all the times)
			cout << "\nQuickHulll\n";
			//change the array of Points to an array of iPairs for use in Quickhull
			iPair a[numberPoints];
			for (int j = 0; j < numberPoints; j++) {
				a[j].first = points[j].x;
				a[j].second = points[j].y;
			}
			//create the clock_t value time
			clock_t time;
			//set the clock
			time = clock();
			//generate the hull using Quickhull
			printHull(a, numberPoints);
			//set time to be the current time minus the start time
			time = clock() - time;
			//output the time in ticks, then seconds
			cout << "completed in " << time << " ticks (" << ((float)time) / CLOCKS_PER_SEC << ")\n";
		}
		return 0;
	}
}
