#include <array>
#include <memory>
#include <iostream>
#include <iterator>
#include <algorithm>
#include <vector>
#include <cassert>
#include <fstream>
#include <string>
#include <chrono>

using namespace std;

int main() {
	std::ifstream myfile("sorted.txt");
	if (myfile.is_open())
	{
		int arrSize = 0;
		int arr[99999]; // Ideally this would be a vector, but you said array

        int x;
	    myfile >> x;
        arr[arrSize++] = x;

		while (true)
		{
            if (myfile.eof())
				break;

            int x;
	        myfile >> x;
            arr[arrSize++] = x;

            if (arr[arrSize - 1] < arr[arrSize - 2]) {
                cout << "Unsorted at " << arrSize << endl;
                return 0;
            }
		}
    // I should have closed the file here, but as the program was ending I was lazy	
    }
	else
	{
		cout << "Unable to open file";
	}
	return 0;
}