#include <iostream>

void Data_recording(string _file_name)
{
	const unsigned int Number_of_Atoms = 10;
	double distance = 2.3;
	string Element_name = "Si";
	double tetta = 2 * M_PI / 5;
	double Radius = distance / 2 / sin(tetta / 2);
	double Height = Radius * cos(tetta / 2);

	double y0 = -Height;
	double sgn = 1;
	double delta_tetta = 0;

	ofstream file;
	file.open(_file_name);
	if (file.is_open())
	{
		file << Number_of_Atoms << endl << endl;
		for (unsigned int i = 0; i < Number_of_Atoms; i++)
		{
			double x = 0, y = y0 - sgn * Radius * cos(tetta * (i + delta_tetta)), z = sgn * Radius * sin(tetta * (i + delta_tetta));
			file.precision(5);
			file << Element_name << "\t\t " << fixed << x << " \t " << y << " \t " << z;

			if (i == 4)
			{
				y0 = Height;
				sgn = -1;
				delta_tetta = -1;
			}
			if (i < Number_of_Atoms - 1) { file << endl; }
		}
		file.close();
	}
	else
	{
		cout << "Error!\t File can't be opened." << endl;
		exit(1);
	}
}
