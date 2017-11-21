#ifndef PARTICLE_H
#define	PARTICLE_H

#define MAX_PARTICLES_IN_ARRAY 2000

enum MISE_pids
{
    kPid_pion = 0,
    kPid_kaon = 1,
    kPid_proton = 2,

    kPid_phi = 3,
    kPid_Lambda = 4,

    kPid_D0 = 5,

//    pid_rho = 3,
};


//particle
struct ParticleDescr
{
    ParticleDescr();
    double eta;
    double phi;
    double pt;
    double ptBeforeKick;
    int charge;
    int pid;
};

//particle array
struct ParticleArr
{
private:
    ParticleDescr fParticles[MAX_PARTICLES_IN_ARRAY];
protected:
    int fNparticles;

public:
    ParticleArr() { fNparticles = 0; }
    int getNparticles() { return fNparticles; }
    //    void setSize( int size ) { fNparticles = size; }
    //    ParticleDescr* getParticle(int id) const { return fNparticles; }

    ParticleDescr *getParticle( int id )
    {
        if ( id >= fNparticles )
            return 0x0;
        return &fParticles[id];
    }
};

#endif	/* PARTICLE_H */

