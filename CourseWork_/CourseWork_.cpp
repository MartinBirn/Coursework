#include <iostream>
#include <ctime>

using namespace std;

enum SOURCE_TYPE
{
	SOURCE_FILE = 1,
	SOURCE_MANUALLY = 2,
	SOURCE_GENERATE = 3,
	SOURCE_EXIT = 4
};

enum FUNC_TYPE
{
	FUNC_DETERMINANT = 1,
	FUNC_INVERSE = 2,
	FUNC_MULTIPLICATION = 3,
	FUNC_EXPONANTION = 4,
	FUNC_STEPPED_APPEARANCE = 5,
	FUNC_EXIT = 6
};

class MatrixClass
{
	int iCurrent_, jCurrent_, iExpected_;
	int rowsCount_, columnsCount_;
	double** matrix_;
	double det_;

public:
	void generateMatrixInMem(int& rowsCount, int& columnsCount, int left_border, int right_border);
	void manualMatrixInput(int& rowsCount, int& columnsCount);
	void readTextMatrix(FILE* stream);
	void countColumns(FILE* stream);
	void countRows(FILE* stream);
	void readNumber(char ch, char* number, int& digit);
	int convert(char* number);
	double computeDet(double** matrix, int rowsCount, double det);
	void findInverseMatrix();
	double** findDownSteppedMatrix(double** matrix, int rowsCount, int columnsCount);
	double** findUpSteppedMatrix(double** matrix, int rowsCount, int columnsCount);
	double* multiDiv(double* str, int columnsCount, double element, double jElement);
	void sumArrs(double*& arr1, double* arr2, int columnsCount, int coeff);
	void addIdentityMatrix(double** matrix);
	void equateDiagonalToUnits(double** matrix, int columnsCount);
	void trimmingTheUnitMatrix(double** matrix, int columnsCount);
	double** matrixMulti(double** matrix1, double** matrix2, int rowsCount1, int columnsCount1, int rowsCount2, int columnsCount2);
	double** exponentMatrix(int exp, int rowsCount, int columnsCount);
	double** getMatrix() { return matrix_; }
	void getParamets(int& rowsCount, int& columnsCount) { rowsCount = rowsCount_; columnsCount = columnsCount_; }
	void setDet(double det) { det_ = det; }
	void memoryClearing(double** matrix, int rowsCount);
};

void outputMatrixToConsole(double** matrix, int rowsCount, int columnsCount)
{
	cout << endl;

	for (int i = 0; i < rowsCount; i++)
	{
		for (int j = 0; j < columnsCount; j++)
		{
			cout << "\t" << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

double MatrixClass::computeDet(double** matrix, int rowsCount, double det)
{
	double** subMatrix;
	int sub_s;

	switch (rowsCount)
	{
	case 1:
		det_ = matrix[0][0];
		return det_;
	case 2:
		det_ = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
		return det_;
	default:
		if (rowsCount < 1)
		{
			cout << "Error." << endl;
			return 0;
		}

		subMatrix = new double* [rowsCount - 1];
		det = 0;

		for (int i = 0; i < rowsCount; i++) //считаем по первому столбцу
		{
			sub_s = 0;

			for (int s = 0; s < rowsCount; s++) //заполняем субматрицу сурс строками
			{
				if (i != s) //исключаем первую строку
				{
					subMatrix[sub_s++] = matrix[s] + 1;
				}
			}

			det = det + pow(-1, i + 2) * matrix[i][0] * MatrixClass::computeDet(subMatrix, rowsCount - 1, det);
		}

		delete[] subMatrix;
		return det;
	}
}

void MatrixClass::findInverseMatrix()
{
	double det_ = computeDet(matrix_, rowsCount_, 0.0);

	if (det_ <= 0.000001 && det_ >= -0.000001)
	{
		cout << "Обратную матрицу найти невозможно, так как определитель равен нулю" << endl;
		exit(2);
	}

	addIdentityMatrix(matrix_); //rowsCount=columnsCount

	matrix_ = findDownSteppedMatrix(matrix_, rowsCount_, columnsCount_ * 2); //columnsCount=rowsCount*2

	matrix_ = findUpSteppedMatrix(matrix_, rowsCount_, columnsCount_ * 2);

	equateDiagonalToUnits(matrix_, columnsCount_ * 2);

	trimmingTheUnitMatrix(matrix_, columnsCount_ * 2);
}

void MatrixClass::addIdentityMatrix(double** matrix)
{
	for (int i = 0; i < rowsCount_; i++)
		for (int j = columnsCount_; j < columnsCount_ * 2; j++)
		{
			if (i + rowsCount_ == j)
				matrix[i][j] = 1;
			else
				matrix[i][j] = 0;
		}

	matrix_ = matrix;

	outputMatrixToConsole(matrix_, rowsCount_, columnsCount_ * 2);
}

double** MatrixClass::findUpSteppedMatrix(double** matrix, int rowsCount, int columnsCount)
{
	for (int i = rowsCount - 2; i > -1; i--)
	{
		for (int j = rowsCount - 1; j > i; j--)
		{
			double* changedArr = multiDiv(matrix[j], columnsCount, matrix[i][j], matrix[j][j]);

			sumArrs(matrix[i], changedArr, columnsCount, -1);
		}
	}

	outputMatrixToConsole(matrix, rowsCount, columnsCount);

	return matrix;
}

double** MatrixClass::findDownSteppedMatrix(double** matrix, int rowsCount, int columnsCount)
{
	for (int i = 1; i < rowsCount; i++)
	{
		for (int j = 0; j < i; j++)
		{
			double* changedArr = multiDiv(matrix[j], columnsCount, matrix[i][j], matrix[j][j]);

			sumArrs(matrix[i], changedArr, columnsCount, -1);
		}
	}

	outputMatrixToConsole(matrix, rowsCount, columnsCount);

	return matrix;
}

double* MatrixClass::multiDiv(double* arr, int columnsCount, double iElement, double jElement)
{
	for (int i = 0; i < columnsCount; i++)
	{
		arr[i] *= iElement / jElement;
	}

	return arr;
}

void MatrixClass::sumArrs(double*& arr1, double* arr2, int columnsCount, int coeff)
{
	for (int k = 0; k < columnsCount; k++)
	{
		arr1[k] += arr2[k] * coeff;

		if (arr1[k]<0.000000001 && arr1[k]>-0.000000001)
		{
			arr1[k] = 0.0;
		}
	}
}

void MatrixClass::equateDiagonalToUnits(double** matrix, int columnsCount)
{
	for (int i = 0; i < rowsCount_; i++)
	{
		double divider = matrix[i][i];

		for (int j = 0; j < columnsCount; j++)
		{
			matrix[i][j] /= divider;
		}
	}

	matrix_ = matrix;

	outputMatrixToConsole(matrix_, rowsCount_, columnsCount);
}

void MatrixClass::trimmingTheUnitMatrix(double** matrix, int columnsCount)
{
	double** tempMatrix = new double* [rowsCount_];

	for (int i = 0; i < rowsCount_; i++) {
		tempMatrix[i] = new double[columnsCount_];
	}

	for (int i = 0, k = 0; i < rowsCount_; i++, k = 0)
	{
		for (int j = columnsCount / 2; j < columnsCount; j++, k++)
		{
			tempMatrix[i][k] = matrix[i][j];
		}
	}

	matrix_ = tempMatrix;
}

double** MatrixClass::matrixMulti(double** matrix1, double** matrix2, int rowsCount1, int columnsCount1, int rowsCount2, int columnsCount2)
{
	double** matrix = new double* [rowsCount1];

	for (int i = 0; i < rowsCount1; i++)
	{
		matrix[i] = new double[columnsCount2];

		for (int j = 0; j < columnsCount2; j++)
		{
			matrix[i][j] = 0;

			for (int k = 0; k < columnsCount1; k++)
				matrix[i][j] += matrix1[i][k] * matrix2[k][j];
		}
	}

	return matrix;
}

double** MatrixClass::exponentMatrix(int exp, int rowsCount, int columnsCount)
{
	double** tempMatrix = matrix_;

	switch (exp)
	{
	case -1:
	{
		findInverseMatrix();

		return matrix_;
		break;
	}
	case 0:
		for (int i = 0; i < rowsCount_; i++)
			for (int j = 0; j < columnsCount_ * 2; j++)
			{
				if (i == j)
				{
					matrix_[i][j] = 1;
				}
				else
				{
					matrix_[i][j] = 0;
				}
			}

		return matrix_;
		break;
	case 1:
		return matrix_;
		break;
	case 2:
		tempMatrix = matrixMulti(matrix_, matrix_, rowsCount, columnsCount, rowsCount, columnsCount);

		return tempMatrix;
		break;
	default:
		if (exp < -1)
		{
			findInverseMatrix();
			tempMatrix = matrix_;

			exp++;

			while (exp++)
			{
				tempMatrix = matrixMulti(tempMatrix, matrix_, rowsCount, columnsCount, rowsCount, columnsCount);
			}
		}
		else if (exp > 2)
		{
			exp--;

			while (exp--)
			{
				tempMatrix = matrixMulti(tempMatrix, matrix_, rowsCount, columnsCount, rowsCount, columnsCount);
			}

			return tempMatrix;
		}
		break;
	}

	return matrix_;
}

void eraseBuffMem()
{
	char c;

	do {
		c = getchar();
	} while (c != '\n' && c != EOF);
}

bool isValidData(char* buff)
{
	int length = strlen(buff);

	for (int i = 0; i < length; i++) {
		if ((buff[i] < 48 || buff[i] > 57) && buff[i] != '-') {
			return false;
		}
	}

	return true;
}

int whatTheSource(int& numberOfSource)
{
	char* buff = new char[10];

	while (true)
	{
		cout << "\nВведите номер источника: ";
		cin >> buff;

		if (isValidData(buff))
		{
			break;
		}
		else
		{
			cout << "Error data. Try again.";
		}
	}

	numberOfSource = atoi(buff);

	delete[] buff;

	return numberOfSource;
}

void MatrixClass::countColumns(FILE* stream) {
	char ch;

	columnsCount_ = 0;

	while (EOF != (ch = fgetc(stream))) {
		if (ch >= 0 || ch <= 9 || ch == '-')
		{
			columnsCount_ = 1;

			break;
		}
	}

	if (!columnsCount_)
	{
		cout << "Error. Matrix has no valid elements" << endl;

		exit(2);
	}

	rewind(stream);

	while (true) {
		ch = fgetc(stream);

		if (ch == ' ') {
			columnsCount_++;
		}
		else if (ch == '\n' || ch == EOF) {
			return;
		}
	}
}

void MatrixClass::countRows(FILE* stream)
{
	rowsCount_ = 1;

	char ch;

	rewind(stream);

	while (EOF != (ch = fgetc(stream))) {
		if (ch == '\n')
		{
			rowsCount_++;
		}
	}
}

int MatrixClass::convert(char* number) {
	int converted = strtol(number, NULL, 10);

	return converted;
}

void MatrixClass::readNumber(char ch, char* number, int& digit)
{
	if ((ch >= '0' && ch <= '9') || ch == '-')
	{
		number[digit++] = ch;
	}
	else if (ch == ' ' || ch == '\n' || ch == EOF)
	{
		if (digit != 0)
		{
			matrix_[iCurrent_][jCurrent_] = MatrixClass::convert(number);

			digit = 0;
		}

		number[0] = '\0';

		jCurrent_++;

		if (ch == '\n' || ch == EOF)
		{
			iCurrent_++;
			jCurrent_ = 0;

			return;
		}

		return;
	}
}

void MatrixClass::readTextMatrix(FILE* stream)
{
	countColumns(stream);
	countRows(stream);

	rewind(stream);

	iCurrent_ = 0, jCurrent_ = 0;

	int digit = 0;

	char number[15];
	char ch;

	matrix_ = new double* [rowsCount_];

	for (int i = 0; i < rowsCount_; i++) {
		matrix_[i] = new double[columnsCount_];
	}

	rewind(stream);

	while (!feof(stream))
	{
		ch = fgetc(stream);

		readNumber(ch, number, digit);
	}
}

void MatrixClass::manualMatrixInput(int& rowsCount, int& columnsCount)
{
	rowsCount_ = rowsCount;
	columnsCount_ = columnsCount;

	matrix_ = new double* [rowsCount_];

	for (int i = 0; i < rowsCount_; i++) {
		matrix_[i] = new double[columnsCount_];
	}

	cout << "Введите элементы матрицы(" << rowsCount << "x" << columnsCount << ")\n";

	for (int i = 0; i < rowsCount; i++)
	{
		for (int j = 0; j < columnsCount; j++)
		{
			cin >> matrix_[i][j];
		}
	}
}

void MatrixClass::generateMatrixInMem(int& rowsCount, int& columnsCount, int left_border, int right_border)
{
	rowsCount_ = rowsCount;
	columnsCount_ = columnsCount;

	matrix_ = new double* [rowsCount_];

	for (int i = 0; i < rowsCount_; i++) {
		matrix_[i] = new double[columnsCount_];
	}

	for (int i = 0; i < rowsCount_; i++) {
		for (int j = 0; j < columnsCount_; j++) {   
			int s;
			s = left_border + rand() % (right_border - left_border + 1);
			matrix_[i][j] = s;
		}
	}

	outputMatrixToConsole(matrix_, rowsCount_, columnsCount_);
}

int whatTheFunc(int& numberOfFunc)
{
	char* buff = new char[10];

	while (true)
	{
		cout << "\nВведите номер действия: ";
		cin >> buff;

		if (isValidData(buff))
		{
			break;
		}
		else
		{
			cout << "Error data. Try again.";
		}
	}

	numberOfFunc = atoi(buff);

	delete[] buff;

	return numberOfFunc;
}

void MatrixClass::memoryClearing(double** matrix, int rowsCount)
{
	for (int i = 0; i < rowsCount; i++)
	{
		delete[] matrix[i];
	}

	delete[] matrix;
}

int main()
{
	srand(time(0));

	setlocale(LC_ALL, "rus");

	int numberOfSource, numberOfFunc;
	int rowsCount, columnsCount, leftBorder, rightBorder;

	double det;

	bool flag;

	MatrixClass obj;

	cout << "Выберите источник матрицы:" << "\n-Чтение матрицы из файла(1)" << "\n-Ввод вручную(2)" << "\n-Генерация матрицы(3)"
		<< "\n-Выйти из программы(4)";

	flag = true;

	while (flag)
	{
		flag = false;

		switch (whatTheSource(numberOfSource))
		{
		case ::SOURCE_FILE:
		{
			FILE* stream;

			char pathToFile[100];

			cout << "Введите путь к файлу, содержащего матрицу, и его название (к примеру C:\matrix.txt)";
			cin >> pathToFile;

			fopen_s(&stream, pathToFile, "r");

			if (stream == NULL)
			{
				perror("Error: ");
				exit(1);
			}

			obj.readTextMatrix(stream);

			double** matrix = obj.getMatrix();
			obj.getParamets(rowsCount, columnsCount);

			outputMatrixToConsole(matrix, rowsCount, columnsCount);

			break;
		}
		case ::SOURCE_MANUALLY:
		{
			cout << "\nВведите размеры матрицы" << "\n-количество строк: ";
			cin >> rowsCount;
			cout << "-количество столбцов: ";
			cin >> columnsCount;

			obj.manualMatrixInput(rowsCount, columnsCount);

			break;
		}
		case ::SOURCE_GENERATE:
			cout << "\nВведите размеры матрицы" << "\n-количество строк: ";
			cin >> rowsCount;
			cout << "-количество столбцов: ";
			cin >> columnsCount;
			cout << "\nВведите границы генерируемых чисел" << "\n-минимальное число:";
			cin >> leftBorder;
			cout << "-максимальное число:";
			cin >> rightBorder;

			obj.generateMatrixInMem(rowsCount, columnsCount, leftBorder, rightBorder);

			break;
		case ::SOURCE_EXIT:
			exit(1);
			break;
		default:
			flag = true;

			cout << "Error number. Try again";

			break;
		}
	}

	while (true)
	{
		cout << "\nВыберите операцию:" << "\n-Найти определитель(1)" << "\n-Найти обратную матрицу(2)" << "\n-Умножение матриц(3)" <<
			"\n-Возвести матрицу в степень(4)" << "\n-Привести матрицу к ступенчатому виду(5)" << "\n-Выйти из программы(6)";

		flag = true;

		while (flag)
		{
			flag = false;

			switch (whatTheFunc(numberOfFunc))
			{
			case ::FUNC_DETERMINANT:
				if (rowsCount != columnsCount)
				{
					flag = true;

					cout << "\nОшибка. Чтобы посчитать определитель, количество строк должно равняться количеству столбцов\n";
				}
				else
				{
					double** matrix = obj.getMatrix();

					det = obj.computeDet(matrix, rowsCount, 0.0);

					obj.setDet(det);

					cout << "\nОпределитель равен = " << endl << det << endl;
				}

				break;
			case ::FUNC_INVERSE:
				if (rowsCount != columnsCount)
				{
					flag = true;

					cout << "\nОшибка. Чтобы найти обратную матрицу, количество строк должно равняться количеству столбцов\n";
				}
				else
				{
					obj.findInverseMatrix();

					double** inverseMatrix = obj.getMatrix();

					outputMatrixToConsole(inverseMatrix, rowsCount, rowsCount);
				}

				break;
			case ::FUNC_MULTIPLICATION:
			{
				MatrixClass obj2;
				int rowsCount2, columnsCount2, leftBorder2, rightBorder2;

				cout << "\nВведите размеры второй матрицы" << "\nколичество строк: ";
				cin >> rowsCount2;
				cout << "количество столбцов: ";
				cin >> columnsCount2;
				cout << "\nВведите границы генерируемых чисел" << "\nминимальное число:";
				cin >> leftBorder2;
				cout << "максимальное число:";
				cin >> rightBorder2;

				obj2.generateMatrixInMem(rowsCount2, columnsCount2, leftBorder2, rightBorder2);
				double** matrix = obj.getMatrix();
				double** matrix2 = obj2.getMatrix();

				if (columnsCount != rowsCount2)
				{
					cout << "Умножение невозможно. Количество столбцов первой матрицы должно ровняться количеству строк второй матрицы" << endl;
				}
				else
				{
					matrix = obj.matrixMulti(matrix, matrix2, rowsCount, columnsCount, rowsCount2, columnsCount2);

					cout << "\nПроизведение матриц равно:\n";
					outputMatrixToConsole(matrix, rowsCount, columnsCount);
				}

				break;
			}
			case ::FUNC_EXPONANTION:
			{
				char* buff = new char[10];
				int exp;

				while (true)
				{
					cout << "\nСтепень возведения матрицы: ";
					cin >> buff;

					if (isValidData(buff))
					{
						break;
					}
					else
					{
						cout << "Error data. Try again.";
					}
				}

				exp = atoi(buff);

				delete[] buff;

				double** matrix = obj.exponentMatrix(exp, rowsCount, columnsCount);

				outputMatrixToConsole(matrix, rowsCount, columnsCount);

				break;
			}
			case ::FUNC_STEPPED_APPEARANCE:
			{
				double** matrix = obj.getMatrix();
				obj.getParamets(rowsCount, columnsCount);

				matrix = obj.findDownSteppedMatrix(matrix, rowsCount, columnsCount);

				outputMatrixToConsole(matrix, rowsCount, columnsCount);

				break;
			}
			case ::FUNC_EXIT:
				exit(1);
				break;
			default:
				flag = true;

				cout << "Error number. Try again";

				break;
			}
		}
	}

	double** matrix = obj.getMatrix();
	obj.getParamets(rowsCount, columnsCount);

	obj.memoryClearing(matrix, rowsCount);

	return 0;
}