// GrapheneAnalyse.cpp : Defines the entry point for the console application.
//

// #include "stdafx.h"
#include "cstdio"
#include <cstring>
#include <cstdlib>
#include <cmath>

struct ATOM {
	int id;
	double z;
	double mass;
	double Ek;
	double Epot;
	int layer;
} *atom;

int main(int argc, char* argv[])
{
	FILE *fp,*fpp;
	double czas, fdummy, z_transm = 15, z_sput = 0., total_organic_mass = 0., total_substrate_mass = 0., total_substrate_sputtered_mass=0., total_projectile_mass = 0.;
	int j, i, id, natoms, nproj = 60, nfile = 0, id_phe = 0, natoms_init = 0, max_layer = 0;
	int n_transm[200][20], n_sput[200],ilayer;
	int n_phe_transm[200], n_phe_sput[200];
	double E_transm[200] , E_sput[200], E_sample[200], Etot[200], Epot[200];
	double E_phe_transm[200], E_phe_sput[200], E_phe_sample[200];
	double d_layer = 3.35;
	int n_transm_proj[200], n_sput_proj[200], n_sample[200], n_phe_sample[200], n_sample_proj[200], n_transm_sum[200], nskip = 0;;
	double E_transm_proj[200], E_sput_proj[200], E_sample_proj[200];
	char hlp[300],elem[4],*p,zbior[300][170];
	int atoms_in_layer = 92162 / 2;
	bool is_q = false;
	atom = nullptr;
	if (argc == 2) {
		sscanf(argv[1], "%d", &nskip);
		printf("%d frames is skipped", nskip);
	}

	system("dir *.lammpstrj > junk");
	fpp = fopen("junk", "r");
	while (!feof(fpp)) {
		fgets(hlp, 300, fpp);
		if (strstr(hlp, "lammpstrj") != nullptr) {
			p = strstr(hlp, "keV");
			strcpy(zbior[nfile], p-2);
			nfile++;
		}
	}
	fclose(fpp);
	for (i = 0; i < 200; i++) {
		for (j = 0; j < 20; j++) n_transm[i][j] = 0;
		n_sput[i]=0;
		n_transm_sum[i] = 0;
		E_transm[i]=E_sput[i]=E_sample[i]=Etot[i]=0.;
		n_phe_transm[i] = n_phe_sput[i] = n_phe_sample[i] =0;
		E_phe_transm[i] = E_phe_sput[i] = E_phe_sample[i] = 0.;
		n_transm_proj[i] = n_sput_proj[i] = n_sample[i] = n_sample_proj[i] = 0;
		E_transm_proj[i] = E_sput_proj[i] = E_sample_proj[i] = Epot[i] = 0.;

	}
	for (j = 0; j < nfile; j++) {
		if ((p = strchr(zbior[j], 10)) != nullptr) *p = 0;
		fp = fopen(zbior[j], "r");
		int len = (int) strlen(zbior[j]);
		if (!fp) {
			printf("File <%s> not found\n", zbior[j]);
			return 1;
		}
		//ITEM: ATOMS id element q x y z mass vx vy vz c_sput_ke c_sput_pe 
		//1 C 4.70215e-06 - 152.018 - 129.85 - 23.9548 12.011 0 0 0 0 - 6.62453
		z_sput = 0;
	e1:	while (!feof(fp)) {
		fgets(hlp, 300, fp);
		if (strstr(hlp, "c_sput_pe") != nullptr) {
			if (strstr(hlp, " q ") != nullptr) is_q = true;
			else is_q = false;
			break;
		}
		if (strstr(hlp, "TIME") != nullptr)
		{
			if (strlen(hlp) > 13) continue;
			fgets(hlp, 300, fp);
			sscanf(hlp, "%lf", &czas);
			continue;
		}
		if (strstr(hlp, "ATOMS") != nullptr) {
			fgets(hlp, 300, fp);
			sscanf(hlp, "%d", &natoms);
			continue;
		}
	}
		if (czas < 1000) {
			int idummy;
			double z;
			if (atom == nullptr) atom = new struct ATOM[natoms+1];
			double Ek_0 = 0., Epot_0 = 0.;
			for (i = 0; i < natoms; i++) {
				fgets(hlp, 300, fp);
				sscanf(hlp, "%d", &id);
				id--;
				if (is_q) sscanf(hlp, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &idummy, elem, &fdummy, &fdummy, &fdummy, &z, &fdummy, &fdummy, &fdummy, &fdummy, &atom[id].Ek, &atom[id].Epot);
				else sscanf(hlp, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf", &idummy, elem, &fdummy, &fdummy, &z, &fdummy, &fdummy, &fdummy, &fdummy, &atom[id].Ek, &atom[id].Epot);
				ilayer = int((-z + d_layer*0.5) / d_layer);
				atom[id].layer = ilayer;
				if (atom[id].Ek == 0 && z_sput > z) z_sput = z;
				Ek_0 += atom[id].Ek;
				Epot_0 += atom[id].Epot;
			}
			goto e1;
		}

		if (nskip > 0) {
			for (i = 0; i < natoms; i++) {
				fgets(hlp, 300, fp);
			}
			nskip--;
			goto e1;
		}
		double z_last = 0;
		for (i = 0; i < natoms; i++) {
			fgets(hlp, 300, fp);
			sscanf(hlp, "%d", &id);
			id--;
			if (is_q) sscanf(hlp, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf", &atom[id].id, elem, &fdummy, &fdummy, &fdummy, &atom[id].z, &atom[id].mass, &fdummy, &fdummy, &fdummy, &atom[id].Ek, &atom[id].Epot);
			else sscanf(hlp, "%d %s %lf %lf %lf %lf %lf %lf %lf %lf %lf", &atom[id].id, elem, &fdummy, &fdummy, &atom[id].z, &atom[id].mass, &fdummy, &fdummy, &fdummy, &atom[id].Ek, &atom[id].Epot);
			Etot[j] += atom[id].Ek;
			Epot[j] += atom[id].Epot;
			if (elem[0] == 'N' && id_phe == 0) id_phe = id;
		}
		fclose(fp);
		natoms_init = natoms;
		natoms -= nproj;
		if (id_phe > 0) natoms = id_phe;
		z_sput -= 10;
		for (i = 0; i < natoms; i++) {
			if (max_layer < atom[i].layer) max_layer = atom[i].layer;
			if (atom[i].z>z_transm) {
				E_transm[j] += atom[i].Ek;
//				if (E_transm[j] > 10000) {
//					printf("zbys %d %.1f %.1f", i, E_transm[j], atom[i].Ek);
//						return 10;
//				}
				ilayer = atom[i].layer;
				n_transm[j][ilayer]++;
				n_transm_sum[j]++;
				total_substrate_mass += atom[i].mass;
			}
			else if (atom[i].z < z_sput) {
				E_sput[j] += atom[i].Ek;
				n_sput[j]++;
				total_substrate_sputtered_mass += atom[i].mass;
			}
			else {
				E_sample[j] += atom[i].Ek;
				n_sample[j]++;
			}
		}
		if (id_phe > 0) { // PHE molecules are present
			for (i = id_phe; i < natoms_init - nproj; i++) {
				if (atom[i].z > z_transm) {
					E_phe_transm[j] += atom[i].Ek;
					n_phe_transm[j]++;
					total_organic_mass += atom[i].mass;
				}
				else if (atom[i].z < z_sput) {
					E_phe_sput[j] += atom[i].Ek;
					n_phe_sput[j]++;
				}
				else {
					E_phe_sample[j] += atom[i].Ek;
					n_phe_sample[j]++;
				}
			}
		}
		for (i = natoms_init-nproj; i < natoms_init; i++) {
			if (atom[i].z>z_transm) {
				E_transm_proj[j] += atom[i].Ek;
				n_transm_proj[j]++;
				total_projectile_mass += atom[i].mass;
			}
			else if (atom[i].z < z_sput) {
				E_sput_proj[j] += atom[i].Ek;
				n_sput_proj[j]++;
			} 
			else {
				E_sample_proj[j] += atom[i].Ek;
				n_sample_proj[j]++;
			}
		}
		printf("File %d\n%d atoms loaded in %d layers\n%d substrate atoms\n%d projectile atoms\n%d organic atoms\nz bottom = %10.1f\n", j + 1, natoms_init, max_layer + 1, natoms, nproj, natoms_init - nproj - id_phe, z_sput);
	}
	printf("%d files analyzed\n\n", nfile);
	double n_transm_ave[20], n_transm_ave_sum = 0 , n_sput_ave = 0;
	double E_transm_ave = 0, E_sput_ave = 0, E_sample_ave = 0, Etot_ave = 0;
	double n_phe_transm_ave = 0, n_phe_sample_ave=0, n_phe_sput_ave = 0;
	double E_phe_transm_ave = 0, E_phe_sput_ave = 0, E_phe_sample_ave = 0, Etot_phe_ave = 0;
	double n_transm_proj_ave = 0, n_sput_proj_ave = 0, n_sample_ave = 0, n_sample_proj_ave = 0;
	double E_transm_proj_ave = 0, E_sput_proj_ave = 0, E_sample_proj_ave = 0;
	for (j = 0; j < nfile; j++) printf("%s %d\n", zbior[j],n_transm[j][0]); printf("\n");
	for (j = 0; j < 20; j++) n_transm_ave[j] = 0;
	for (i = 0; i < nfile; i++) {
		E_transm_ave += E_transm[i];
		E_sput_ave += E_sput[i]; 
		E_sample_ave += E_sample[i];
		Etot_ave += Etot[i];
		for (j = 0; j < 20; j++) n_transm_ave[j] += n_transm[i][j];
		n_sput_ave += n_sput[i];
		n_sample_ave += n_sample[i];

		E_phe_transm_ave += E_phe_transm[i];
		E_phe_sput_ave += E_phe_sput[i];
		E_phe_sample_ave += E_phe_sample[i];
//		Etot_phe_ave += Etot_phe[i];
		n_phe_transm_ave += n_phe_transm[i];
		n_phe_sput_ave += n_phe_sput[i];
		n_phe_sample_ave += n_phe_sample[i];

		E_transm_proj_ave += E_transm_proj[i];
		E_sput_proj_ave += E_sput_proj[i];
		E_sample_proj_ave += E_sample_proj[i];
		n_transm_proj_ave += n_transm_proj[i];
		n_sput_proj_ave += n_sput_proj[i];
		n_sample_proj_ave += n_sample_proj[i];
	}
	for (j = 0; j < 20; j++) n_transm_ave_sum += n_transm_ave[j];
	E_transm_ave /= double(nfile);
	E_sput_ave /= double(nfile);
	E_sample_ave /= double(nfile);
	Etot_ave /= double(nfile);
	for (j = 0; j < 20; j++) n_transm_ave[j] /= double(nfile);
	n_transm_ave_sum /= double(nfile);
	n_sput_ave /= double(nfile);
	n_sample_ave /= double(nfile);

	E_phe_transm_ave /= double(nfile);
	E_phe_sput_ave /= double(nfile);
	E_phe_sample_ave /= double(nfile);
//	Etot_ave /= double(nfile);
	n_phe_transm_ave /= double(nfile);
	n_phe_sput_ave /= double(nfile);
	n_phe_sample_ave /= double(nfile);

	E_transm_proj_ave /= double(nfile);
	E_sput_proj_ave /= double(nfile);
	E_sample_proj_ave /= double(nfile);
	n_transm_proj_ave /= double(nfile);
	n_sput_proj_ave /= double(nfile);
	n_sample_proj_ave /= double(nfile);
	total_organic_mass /= double(nfile);
	total_substrate_mass /= double(nfile);
	total_projectile_mass /= double(nfile);
	total_substrate_sputtered_mass /= double(nfile);

	double sn_transm_ave[20], sn_transm_ave_sum = 0, sn_sput_ave = 0;
	double sE_transm_ave = 0, sE_sput_ave = 0, sE_sample_ave = 0, sEtot_ave = 0;
	double sn_phe_transm_ave = 0, sn_phe_sput_ave = 0, sn_phe_sample_ave = 0;
	double sE_phe_transm_ave = 0, sE_phe_sput_ave = 0, sE_phe_sample_ave = 0;
	double sn_transm_proj_ave = 0, sn_sput_proj_ave = 0, sn_sample_ave = 0, sn_sample_proj_ave = 0;
	double sE_transm_proj_ave = 0, sE_sput_proj_ave = 0, sE_sample_proj_ave = 0;
	for (j = 0; j < 20; j++) sn_transm_ave[j] = 0;
	for (i = 0; i < nfile; i++) {
		sE_transm_ave += (E_transm[i]- E_transm_ave)*(E_transm[i] - E_transm_ave);
		sE_sput_ave += (E_sput[i]- E_sput_ave)*(E_sput[i] - E_sput_ave);
		sE_sample_ave += (E_sample[i]- E_sample_ave)*(E_sample[i] - E_sample_ave);
		sEtot_ave += (Etot[i]- Etot_ave)*(Etot[i] - Etot_ave);
		for (j = 0; j < 20; j++) {
			sn_transm_ave[j] += (n_transm_ave[j] - n_transm[i][j])*(n_transm_ave[j] - n_transm[i][j]);
		}
		sn_transm_ave_sum += (n_transm_ave_sum - n_transm_sum[i])*(n_transm_ave_sum - n_transm_sum[i]);
		sn_sput_ave += (n_sput[i]- n_sput_ave)*(n_sput[i] - n_sput_ave);
		sn_sample_ave += (n_sample[i]- n_sample_ave)*(n_sample[i] - n_sample_ave);

		sE_phe_transm_ave += (E_phe_transm[i] - E_phe_transm_ave)*(E_phe_transm[i] - E_phe_transm_ave);
		sE_phe_sput_ave += (E_phe_sput[i] - E_phe_sput_ave)*(E_phe_sput[i] - E_phe_sput_ave);
		sE_phe_sample_ave += (E_phe_sample[i] - E_phe_sample_ave)*(E_phe_sample[i] - E_phe_sample_ave);
		sn_phe_transm_ave += (n_phe_transm[i] - n_phe_transm_ave)*(n_phe_transm[i] - n_phe_transm_ave);
		sn_phe_sput_ave += (n_phe_sput[i] - n_phe_sput_ave)*(n_phe_sput[i] - n_phe_sput_ave);
		sn_phe_sample_ave += (n_phe_sample[i] - n_phe_sample_ave)*(n_phe_sample[i] - n_phe_sample_ave);

		sE_transm_proj_ave += (E_transm_proj[i]- E_transm_proj_ave)*(E_transm_proj[i] - E_transm_proj_ave);
		sE_sput_proj_ave += (E_sput_proj[i]- E_sput_proj_ave)*(E_sput_proj[i] - E_sput_proj_ave);
		sE_sample_proj_ave += (E_sample_proj[i]- E_sample_proj_ave)*(E_sample_proj[i] - E_sample_proj_ave);
		sn_transm_proj_ave += (n_transm_proj[i]- n_transm_proj_ave)*(n_transm_proj[i] - n_transm_proj_ave);
		sn_sput_proj_ave += (n_sput_proj[i]- n_sput_proj_ave)*(n_sput_proj[i] - n_sput_proj_ave);
		sn_sample_proj_ave += (n_sample_proj[i]- n_sample_proj_ave)*(n_sample_proj[i] - n_sample_proj_ave);
	}

	sE_transm_ave /= double(nfile)*double(nfile - 1);
	sE_sput_ave /= double(nfile)*double(nfile - 1);
	sE_sample_ave /= double(nfile)*double(nfile - 1);
	sEtot_ave /= double(nfile)*double(nfile - 1);
	for (j = 0; j < 20;j++) sn_transm_ave[j] /= double(nfile)*double(nfile - 1);
	sn_sput_ave /= double(nfile)*double(nfile - 1);
	sn_sample_ave /= double(nfile)*double(nfile - 1);
	sn_transm_ave_sum /= double(nfile)*double(nfile - 1);

	sE_phe_transm_ave /= double(nfile)*double(nfile - 1);
	sE_phe_sput_ave /= double(nfile)*double(nfile - 1);
	sE_phe_sample_ave /= double(nfile)*double(nfile - 1);
	sn_phe_transm_ave /= double(nfile)*double(nfile - 1);
	sn_phe_sput_ave /= double(nfile)*double(nfile - 1);
	sn_phe_sample_ave /= double(nfile)*double(nfile - 1);

	sE_transm_proj_ave /= double(nfile)*double(nfile - 1);
	sE_sput_proj_ave /= double(nfile)*double(nfile - 1);
	sE_sample_proj_ave /= double(nfile)*double(nfile - 1);
	sn_transm_proj_ave /= double(nfile)*double(nfile - 1);
	sn_sput_proj_ave /= double(nfile)*double(nfile - 1);
	sn_sample_proj_ave /= double(nfile)*double(nfile-1);

	sE_transm_ave = sqrt(sE_transm_ave);
	sE_sput_ave = sqrt(sE_sput_ave);
	sE_sample_ave = sqrt(sE_sample_ave);
	sEtot_ave = sqrt(sEtot_ave);
	for (j = 0; j < 20; j++) sn_transm_ave[j] = sqrt(sn_transm_ave[j]);
	sn_transm_ave_sum = sqrt(sn_transm_ave_sum);
	sn_sput_ave = sqrt(sn_sput_ave);
	sn_sample_ave = sqrt(sn_sample_ave);

	sE_phe_transm_ave = sqrt(sE_transm_ave);
	sE_phe_sput_ave = sqrt(sE_sput_ave);
	sE_phe_sample_ave = sqrt(sE_sample_ave);
	sn_phe_transm_ave = sqrt(sn_phe_transm_ave);
	sn_phe_sput_ave = sqrt(sn_sput_ave);
	sn_phe_sample_ave = sqrt(sn_sample_ave);

	sE_transm_proj_ave = sqrt(sE_transm_proj_ave);
	sE_sput_proj_ave = sqrt(sE_sput_proj_ave);
	sE_sample_proj_ave = sqrt(sE_sample_proj_ave);
	sn_transm_proj_ave = sqrt(sn_transm_proj_ave);
	sn_sput_proj_ave = sqrt(sn_sput_proj_ave);
	sn_sample_proj_ave = sqrt(sn_sample_proj_ave);
	 
	fclose(fpp);
	printf("\nSample atoms\n\nForward ejected = %.1f +- %.1f\nSputtered = %.1f +- %.1f\nIn Sample = %.1f +- %.1f\nForward ejected substrate mass = %.1f\nSputtered substrate mass = %.1f\n", n_transm_ave_sum, sn_transm_ave_sum, n_sput_ave, sn_sput_ave, n_sample_ave, sn_sample_ave,total_substrate_mass, total_substrate_sputtered_mass);
	printf("\nKE          \nForward ejected = %.1f +- %.1f\nSputtered = %.1f +- %.1f\nIn Sample = %.1f +- %.1f\n", E_transm_ave, sE_transm_ave, E_sput_ave, sE_sput_ave, E_sample_ave, sE_sample_ave);
	printf("\n\nProjectile atoms\n\nTransmitted = %.1f +- %.1f\nBackreflected = %.1f +- %.1f\nIn Sample = %.1f +- %.1f\nEjected mass = %.1f\n", n_transm_proj_ave, sn_transm_proj_ave, n_sput_proj_ave, sn_sput_proj_ave, n_sample_proj_ave, sn_sample_proj_ave,total_projectile_mass);
	printf("\nKE          \nTransmitted = %.1f +- %.1f\nBackreflected = %.1f +- %.1f\nIn Sample = %.1f +- %.1f\n", E_transm_proj_ave, sE_transm_proj_ave, E_sput_proj_ave, sE_sput_proj_ave, E_sample_proj_ave, sE_sample_proj_ave);
	if (id_phe > 0) {
		printf("\nOrganic atoms\n\nForward ejected = %.1f +- %.1f\nSputtered = %.1f +- %.1f\nIn Sample = %.1f +- %.1f\nEjected mass = %.1f\n", n_phe_transm_ave, sn_phe_transm_ave, n_phe_sput_ave, sn_phe_sput_ave, n_phe_sample_ave, sn_phe_sample_ave, total_organic_mass);
		printf("\nKE           \nForward ejected = %.1f +- %.1f\nSputtered = %.1f +- %.1f\nIn Sample = %.1f +- %.1f\n", E_phe_transm_ave, sE_phe_transm_ave, E_phe_sput_ave, sE_phe_sput_ave, E_phe_sample_ave, sE_phe_sample_ave);
	}
	printf("\n\nEtot = %.1f +- %1f\nTotal ejected mass = %.1f\n", Etot_ave, sEtot_ave, total_organic_mass + total_projectile_mass + total_substrate_mass);

	p = strchr(zbior[0], 'L'); p++; *p = 0; strcat(zbior[0], ".dat");
	fp = fopen(zbior[0], "w");
	fprintf(fp,"%d atoms loaded\n%d projectile atoms\nz bottom = %10.1f\n", natoms, nproj, z_sput);
	fprintf(fp,"%d files analyzed\n\n", nfile);
	fprintf(fp, "\nSample atoms\n\nForward ejected = %.1f +- %.1f\nSputtered = %.2f +- %.2f\nIn Sample = %.1f +- %.1f\nForward ejected substrate mass = %.1f\nSputtered substrate mass = %.1f\n", n_transm_ave_sum, sn_transm_ave_sum, n_sput_ave, sn_sput_ave, n_sample_ave, sn_sample_ave,total_substrate_mass, total_substrate_sputtered_mass);
	fprintf(fp, "\nKE          \nForward ejected = %.1f +- %.1f\nSputtered = %.2f +- %.2f\nIn Sample = %.1f +- %.1f\n", E_transm_ave, sE_transm_ave, E_sput_ave, sE_sput_ave, E_sample_ave, sE_sample_ave);
	fprintf(fp, "\n\nProjectile atoms\n\nTransmitted = %.1f +- %.1f\nBackreflected = %.2f +- %.2f\nIn Sample = %.1f +- %.1f\nEjected mass = %.1f\n", n_transm_proj_ave, sn_transm_proj_ave, n_sput_proj_ave, sn_sput_proj_ave, n_sample_proj_ave, sn_sample_proj_ave,total_projectile_mass);
	fprintf(fp, "\nKE          \nTransmitted = %.1f +- %.1f\nBackreflected = %.2f +- %.2f\nIn Sample = %.1f +- %.1f\n", E_transm_proj_ave, sE_transm_proj_ave, E_sput_proj_ave, sE_sput_proj_ave, E_sample_proj_ave, sE_sample_proj_ave);
	if (id_phe > 0) {
		fprintf(fp, "\nOrganic atoms\n\nForward ejected = %.1f +- %.1f\nSputtered = %.2f +- %.2f\nIn Sample = %.1f +- %.1f\nEjected_mass = %.1f\n", n_phe_transm_ave, sn_phe_transm_ave, n_phe_sput_ave, sn_phe_sput_ave, n_phe_sample_ave, sn_phe_sample_ave, total_organic_mass);
		fprintf(fp,"\nKE           \nForward ejected = %.1f +- %.1f\nSputtered = %.2f +- %.2f\nIn Sample = %.1f +- %.1f\n", E_phe_transm_ave, sE_phe_transm_ave, E_phe_sput_ave, sE_phe_sput_ave, E_phe_sample_ave, sE_phe_sample_ave);
	}
	fprintf(fp, "\n\nEtot = %.1f +- %1f\nTotal ejected mass = %.1f\n", Etot_ave, sEtot_ave, total_organic_mass + total_projectile_mass + total_substrate_mass);
	fprintf(fp, "\n");
	for (i = 0; i <= max_layer; i++) fprintf(fp, "%d %10.3f +/- %10.3f\n", i + 1, n_transm_ave[i],sn_transm_ave[i]);
	fprintf(fp, "\n");
	for (i = 0; i < nfile; i++) fprintf(fp, "%d %d\n", i + 1, n_transm_sum[i]);  
	fclose(fp);
	return 0;
}

