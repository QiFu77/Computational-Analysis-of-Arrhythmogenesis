
/**
 * Construct geometry from geometry files
 *
 */
int*** constructGeometry()
{
	int x,y,z;// cell loop index
	// Stupid C++ can only initialize a multi-dimen array like this. If not, the array cannot be returned.
	int ***geo = new int** [NX];
	for (x = 0; x < NX; x++)
	{
		geo[x] = new int*[NY];
		for (y = 0; y < NY; y++)
		{
			geo[x][y] = new int[NZ];
			for (z = 0; z < NZ; z++)
				geo[x][y][z] = 0;
		}
	}

	printf("Construct geometry from files...\n");
	FILE *geoFile = fopen("Geometry/BensonHumanVentricle_102*102*102/heterogeneity.txt", "rt");
	
	
	int c;// cell type
	int count = 0;
	for (z = 0; z < NZ; z++) 
	{// be careful about the reading sequence
		for (y = 0; y < NY; y++) 
		{
			for (x = 0; x < NX; x++) 
			{ 
				fscanf(geoFile, "%d ", &c);
				if (c == 0) geo[x][y][z] = 0;	// node is outside tissue
				if (c == 1) geo[x][y][z] = 1;   // node is endocardial surface (for stimulation)
				if (c == 2) geo[x][y][z] = 2;	// node is endocardial
				if (c == 3) geo[x][y][z] = 3;	// node is midmyocardial
				if (c == 4) geo[x][y][z] = 4;	// node is epicardial

				if(x == 0 || x== NX-1 || y == 0 || y==NY-1 || z == 0 || z == NZ -1)
					if (geo[x][y][z] != 0)  cout << "Tissue exists in geometry boundary!" << endl;
				if (c != 0) count++;
			}
			fscanf(geoFile, "\n");
		}
		fscanf(geoFile, "\n");
	}
	fclose(geoFile);
	cout << "Tissue cell number = " << count << endl;
	return geo;
}
