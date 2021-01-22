#include <iostream>
#include <cstring>
#include<fstream>
#include<math.h>
#define MAX 1000
using namespace std;

int GetMax(int N1,int N2,int N3); // function to get max from three numbers
float CalculateAlignment(string sequence1,string sequence2);		// function to set score matrix
float CalculatePercentage(string seq1, string aligned);		

template <class T>
class DynamicSafeArray2D
{
	private : 
				T **data;
				int rows,cols;
	public : 
			DynamicSafeArray2D()
			{
				data = 0 ; 
				rows = 0 ; 
				cols = 0 ; 	
			}
			DynamicSafeArray2D(int r,int c)
			{
				rows = r;
				cols = c;
				data = new T *[rows];
					
				for(int i = 0; i < rows; i++)
				{
					data[i] = new T[cols] ;				
				}	
			}
			T& operator ()(int r,int c)
			{
				if(r>=0 && r<= rows && c>=0 && c<=cols)
				{
					return data[r][c];	
				}	
			}
			void display()
			{
				for(int i=0;i<rows;i++)
				{
					for(int j=0;j<cols;j++)
					{	
						cout<<data[i][j]<<"\t";
					}
					cout<<endl;		
				}	
			}
			DynamicSafeArray2D (const DynamicSafeArray2D &rhs)
			{
				rows = rhs.rows;
				cols = rhs.cols;
				data = new T*[rows];	
				int i,j ;
				for(i = 0 ; i < rows ; i++ )
				{
					data[i] = new T[cols];
				}
				for(i = 0; i < rows; i++)
				{																
					for(j = 0; j < cols; j++)				
					{
						data[i][j] = rhs.data[i][j];						// copying data 
					}
				}
			}
			DynamicSafeArray2D& operator = (const DynamicSafeArray2D& rhs)
			{
				if (this == &rhs)
       			{
      				return *this;		
				}	
				for (int i = rows-1; i >=0 ; i--)
				{
					delete data[i];
				}
				delete [] data;
	
				rows = rhs.rows;
				cols = rhs.cols;
				data = new T*[rows];
				
				int i, j ; 
				for(i = 0 ; i < rows ; i++ )
				{
					data[i] = new T[cols];
				}
				
				for(i = 0; i < rows; i++)
				{
					for(j = 0; j < cols; j++)
					{
						data[i][j] = rhs.data[i][j];	// copying data 
					}
				}				
				return *this;
			}		
};
int main()
{	
	float p1,p2,p3,max;
	ifstream in;
	ofstream out;
	
	string S1 ;
	string S2;
	cin>>S2;
	in.open("file1.txt");			// opeing file 1
	in>>S1;							// reading data base sequence 1
	in.close();						// closing file
	p1=CalculateAlignment(S2,S1);
	cout<<endl;
	
	
	in.open("file2.txt");			// opeing file 2
	in>>S1;							// reading data base sequence 2
	in.close();						// closing file	
	p2=CalculateAlignment(S2,S1);
	cout<<endl;
	
	in.open("file3.txt");			// opeing file 3
	in>>S1;							// reading data base sequence 3
	in.close();						// closing file	
	p3=CalculateAlignment(S2,S1);
	cout<<endl;
	
	
	if(p1>p2 && p1>p3)
	{
		cout<<"MyoGlobin"<<endl;
	}
	else if(p2>p1 && p2>p3)
	{
		cout<<"Lysozyme"<<endl;
	}
	else if(p3>p1 && p3>p2)
	{
		cout<<"Protein-Ribonuclease-Seminal"<<endl;
	}
	
}
int GetMax(int N1,int N2,int N3)
{
	int Max;
	if(N1>=N2 && N1>=N3)
	{
		Max = N1 ;
	}
	else if(N2>=N1 && N2>=N3)
	{
		Max = N2 ;
	}
	else if(N3>=N1 && N3>=N2)
	{
		Max = N3 ;
	}
	return Max;
}


float CalculateAlignment(string sequence1,string sequence2)
{
	int gape = -1;				// for a gape we add -1 as penality
	int mismatch = -1;			// for a mismatch we add -1 as penality
	int match = 1 ;				// for a match we add 1 to score 	
	
	//int cnt=0;
	int LengthOfS1 = sequence1.length()+1;			// lengthof sequence 1
	int LengthOfS2 = sequence2.length()+1;			// lengthof sequence 2
	
	cout<<"\n\nOriginal Sequences :\n";
	cout<<sequence2<<endl;
	cout<<sequence1<<endl;
	
	DynamicSafeArray2D<int> ScoreMatrix(LengthOfS1,LengthOfS2);
	DynamicSafeArray2D<int> TraceMatrix(LengthOfS1,LengthOfS2);

	
	int i , j ,top=0,left=0,diag=0;
	for( i = 0 ; i < LengthOfS2; i++)
	{
		ScoreMatrix(0,i)= gape*i;// adding gape to first row 
		TraceMatrix(0,i)=1;
	}
	for( i = 0 ; i < LengthOfS1; i++)
	{
		ScoreMatrix(i,0) = gape*i;	// adding gape to first column 
		TraceMatrix(i,0)=-1;
	}
	TraceMatrix(0,0)=0;
	for( i = 1; i < LengthOfS1; i++)
	{
		for(j = 1; j < LengthOfS2; j++)
		{						
			top = ScoreMatrix(i-1,j)+gape;							// Top element of i,j position in 2D Matrix 
			left = ScoreMatrix(i,j-1)+gape;							// left element of i,j position in 2D Matrix
			diag = ScoreMatrix(i-1,j-1);							// Top Left (Diagonal) element of i,j position in 2D Matrix
			
			if(sequence1[i-1] == sequence2[j-1])					
			{
				diag = diag+match;
			}
			else
			{
				diag = diag+mismatch;		// since gape and mistach have same score that is -1 
			}
			
			int temp = GetMax(top,left,diag); 
			if( temp == top )	//up so -1
			{			
				ScoreMatrix(i,j) = top ;
				TraceMatrix(i,j)=-1;		
			}
			else if( temp == left)	//left so 1
			{
				ScoreMatrix(i,j) = left ;
				TraceMatrix(i,j)=1;	
			}
			else if( temp == diag)	//diagonal so 0
			{
				ScoreMatrix(i,j) = diag ;
				TraceMatrix(i,j)=0;	
			}			
		}
	}
	
	ScoreMatrix.display();
	cout<<endl;
	
	cout<<"Aligned Sequences:\n";
	cout<<sequence2<<endl;		// printing database sequence
	
	i = LengthOfS1-1;
	j = LengthOfS2-1;
	string alignedsequence ="";
	// tracing back
	while(i!=0 && j!=0)
	{
			top = ScoreMatrix(i-1,j);							// Top element of i,j position in 2D Matrix 
			left = ScoreMatrix(i,j-1);							// left element of i,j position in 2D Matrix
			diag = ScoreMatrix(i-1,j-1);
			
			if(TraceMatrix(i,j)==0)	
			{
				alignedsequence = alignedsequence+sequence1[i-1];
				i = i - 1;
			 	j = j - 1;
			 	//cnt++;
			}
			else if(TraceMatrix(i,j)==-1)	
			{
				alignedsequence = alignedsequence+sequence1[i-1];
				i = i - 1;
				//cnt++;
			}
			else
			{
				alignedsequence = alignedsequence+"-";
				j--;
				
			}
		
	}
	string alignedsequence2="" ;						
	for(int i = alignedsequence.length()-1 ;i>=0 ;i--)
	{
		alignedsequence2 =alignedsequence2+ alignedsequence[i] ;			// reversing

	}
	cout<<alignedsequence2;	
	
	float ret=CalculatePercentage(sequence1, alignedsequence2);
	return ret;
}

float CalculatePercentage(string seq1, string aligned)
{
	float count=0;
	int i,j;
	int LengthOfS1 = seq1.length()+1;		
	int LengthOfAligned = aligned.length()+1;		

	for( i = 0; i < LengthOfS1; i++)
	{
		if(seq1[i] == aligned[i])					//(matchedCount/LengthOfS1)*100
		{
			count++;
		}		
	}
	cout<<endl<<count<<endl;
	float ans=(count/LengthOfS1) *100;
	
	cout<<endl<<"Percentage: "<<ans;
	
	return ans;
}
