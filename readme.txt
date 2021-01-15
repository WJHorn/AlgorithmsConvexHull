this is the readme file for Will Horn's Algorithms project 2 (Convex Hull)

main.cpp:
	main.cpp includes implementations for three different ConvexHull Algorithms
	1) Grahm Scan:
		grahm scan finds the bottom-most point (if there are two bottom-most points, then use the bottom-most, left-most point) and makes it p0 (first point in the hull), then
		it considers the remaining n-1 points and sorts them by angle (counterclockwise) in relation to p0
		after it sorts them, after sorting, it checks if two or more points have the same angle, and removes all but the furthest point with the same angle.
		next it creates an empty stack and pushes points p0, p1, and p2 onto it
		then it processes the remaining points one by one for every point
		if it doesn't make a left turn, it removes the previous 2 points from the stack
		lastly, once all hull points have been pushed onto the stack, it saves all of these points to the file hull_G.txt

	2) Jarvis March:
		Jarvis march initializes p as the leftmost point
		the rest is done while we have not returned to the first point
			find the point that makes the smallest angle with y = p and p->p1 where p1 is the afformentioned point
			store p1 as the next point in the convex hull
			find the point that makes the smallest angle with p->p1 and p1->p2
			store p2 as the next point in the convex hull
			...
		we have reached p again
		Jarvis march now stores the values of the hull to the file hull_J.txt

	3) Quickhull:
		quickhull finds the point with the minimum x-coord and the maximum x-coord
		it then makes a line joining the two points
		next, it finds the point with the maximum distance from this line, and form a triangle with the max-x and min-x
			anything in this triangle will not be part of the convex hull
		the previus step divides into 2 sub-problems that are recursively called in order to find all of the points on the hull
		once there are no more points on the outside of the hull, then the values are inserted into the hull
		finally, quickhull saves the hull values to the file hull_Q.txt

	RUNNING MAIN.CPP:
		once compiled using the linux command "g++ main.cpp", you can run the program using one of three options
		(these are based on running them on the linux terminal for the lab machines)
			"./a.out G test.txt"
				runs using Grahm Scan
			"./a.out J test.txt"
				runs using Jarvis March
			"./a.out Q test.txt"
				runs using Quickhull

	ISSUES WITH MAIN.CPP:
		There was only one problem that I ran into. This was getting the hull file to be in the correct order for quicksort
		The implementation I used used a set to hold all of the points, and thereby it sorted all of them, but also removed duplicates.
		It runs in the correct time, which I think is the most impoirtant part, but the hull in the Gui goes in numerical order, so it is not fully correct.
		I attempted to use a vector<Point> rather than a set of points, but it had a lot of duplicates, which made the hull have lines going along the lines 
		that were created in previous steps.