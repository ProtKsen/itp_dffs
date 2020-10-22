#include "Tests.h"

#include "CVs/RMSDCV.h"
#include <iostream>
#include <fstream>

using namespace SSAGES;

class RMSDCVTests : public ::testing::Test {
protected:
	virtual void SetUp()
	{
		//Set up xyz files
		std::ofstream myfile;

		//Normal file with extra line (should pass)
		myfile.open (filexyz1);
		myfile << "3\n";
		myfile << "Comments here\n";
		myfile << "1 0.0 2.2 1.2\n";
		myfile << "2 1.2 2.5 0.5\n";
		myfile << "3 1.2 1.0 1.3";
		myfile.close();

		normal_RMSD = new RMSDCV({1,3}, filexyz1, true);
		no_range_RMSD = new RMSDCV({1,2,3}, filexyz1, false);

		// Set up snapshot No. 1
		// Snapshot 1 contains three atoms
		snapshot1 = new Snapshot(comm, 0);

		Matrix3 H;
		H << 100.0, 0.0, 0.0,
		     0.0, 100.0, 0.0,
		     0.0, 0.0, 100.0;
		snapshot1->SetHMatrix(H);

		unsigned int n = 3;
		snapshot1->SetNumAtoms(n);

		auto& pos1 = snapshot1->GetPositions();
		pos1.resize(n);
		pos1[0][0] = 0.0;
		pos1[0][1] = 2.2;
		pos1[0][2] = 0.9;

		pos1[1][0] = 1.2;
		pos1[1][1] = 2.1;
		pos1[1][2] = 0.4;

		pos1[2][0] = 1.2;
		pos1[2][1] = 1.1;
		pos1[2][2] = 1.1;

		auto& ids1 = snapshot1->GetAtomIDs();
		ids1.resize(n);
		ids1[0] = 1;
		ids1[1] = 2;
		ids1[2] = 3;

		auto& mass = snapshot1->GetMasses();
		mass.resize(n);
		mass[0] = 1;
		mass[1] = 1;
		mass[2] = 1;
	}

	virtual void TearDown()
	{
		std::remove(filexyz1.c_str());
		delete normal_RMSD;
		delete no_range_RMSD;
		delete snapshot1;
	}

	std::string filexyz1 = "test1.xyz";
	RMSDCV *normal_RMSD;
	RMSDCV *no_range_RMSD;

	mxx::comm comm;

	Snapshot *snapshot1;
};

TEST_F(RMSDCVTests, RMSDTest)
{
	normal_RMSD->Initialize(*snapshot1);
	normal_RMSD->Evaluate(*snapshot1);
	EXPECT_NEAR(normal_RMSD->GetValue(), 0.20188203060665863, 1000*eps);
}

TEST_F(RMSDCVTests, RMSDRangeTest)
{
	normal_RMSD->Initialize(*snapshot1);
	normal_RMSD->Evaluate(*snapshot1);
	no_range_RMSD->Initialize(*snapshot1);
	no_range_RMSD->Evaluate(*snapshot1);

	EXPECT_EQ(normal_RMSD->GetValue(), no_range_RMSD->GetValue());
}

TEST_F(RMSDCVTests, RMSDBadRangeTest)
{
	EXPECT_ANY_THROW(new RMSDCV({1,2,3}, filexyz1, true));
	EXPECT_ANY_THROW(new RMSDCV({3,1}, filexyz1, true));
}
