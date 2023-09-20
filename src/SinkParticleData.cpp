#include <AMRLevelOrion.H>
#include <SinkParticleData.H>

//#include "ParticleList.H"

#define IOPREC 20

int SinkParticleData::angmom_method = 0;

ostream &operator<<(ostream &os, const SinkParticleData& partdata) {
  // Set precision to IOPREC decimal places
  int oldprec = os.precision(IOPREC);

  // Output order is mass, position, momentum
  os << partdata.m << " ";
  for (int i=0; i<CH_SPACEDIM; i++)
    os << partdata.pos[i] << " ";
  for (int i=0; i<CH_SPACEDIM; i++)
    os << partdata.mom[i] << " ";
  for (int i=0; i<CH_SPACEDIM; i++)
    os << partdata.angmom[i] << " ";

  // Restore old precision
  os.precision(oldprec);

  return(os);
}

istream &operator>>(istream &is, SinkParticleData& partdata) {
  // Set precision to IOPREC decimal places
  int oldprec = is.precision(IOPREC);

  // Input order is mass, position, momentum
  is >> partdata.m;
  for (int i=0; i<CH_SPACEDIM; i++) is >> partdata.pos[i];
  for (int i=0; i<CH_SPACEDIM; i++) is >> partdata.mom[i];
  for (int i=0; i<CH_SPACEDIM; i++) is >> partdata.angmom[i];

  // Restore old precision
  is.precision(oldprec);

  // Particles that have been read in are not new
#ifdef MERGENEW
  partdata.newParticle = 0;
#endif

  return(is);
}

SinkParticleData::SinkParticleData(Real *newpos, Real *newmom,
				   Real *newangmom, Real newmass) :
  ParticleData(newpos)
{
  // get the curent angular momentum treatment option
  for (int i=0; i<CH_SPACEDIM; i++) mom[i]=newmom[i];
  for (int i=0; i<CH_SPACEDIM; i++) angmom[i]=newangmom[i];
  m = newmass;
  m_wind_merge=0;
#ifdef MERGENEW
  // This particle is new
  newParticle=1;
#endif
}

SinkParticleData::~SinkParticleData() {}

SinkParticleData&
SinkParticleData::operator=(const SinkParticleData& rhs) {
  for (int i=0; i<CH_SPACEDIM; i++) {
    pos[i]=rhs.pos[i];
    mom[i]=rhs.mom[i];
    angmom[i]=rhs.angmom[i];
  }
  m = rhs.m;
#ifdef MERGENEW
  newParticle = rhs.newParticle;
#endif
  return(*this);
}

SinkParticleData
SinkParticleData::operator+(const SinkParticleData& rhs) {
  SinkParticleData result;
  result.m = m+rhs.m;
  for (int i=0; i<CH_SPACEDIM; i++) {
    result.pos[i]=(pos[i]*m+rhs.pos[i]*rhs.m)/result.m;
    result.mom[i]=mom[i]+rhs.mom[i];
    result.angmom[i]=angmom[i]+rhs.angmom[i];
  }
  // Add the orbital angular momentum of the two stars about the center
  // of mass to the spin angular momentum of the individual stars
  if (CH_SPACEDIM==3) {
    Real r1[CH_SPACEDIM], r2[CH_SPACEDIM], deltal1[CH_SPACEDIM], deltal2[CH_SPACEDIM], p1[CH_SPACEDIM], p2[CH_SPACEDIM];
    for (int i=0; i<CH_SPACEDIM; i++) {
      deltal1[i] = 0.;
      deltal2[i] = 0.;
      r1[i] = pos[i] - result.pos[i];
      r2[i] = rhs.pos[i] - result.pos[i];
      // momentum in mesh frame
      p1[i]=mom[i];
      p2[i]=rhs.mom[i];
      // put the momentum in the frame of the centroid
      p1[i]=m*(p1[i]/m - result.mom[i]/result.m);
      p2[i]=rhs.m*(p2[i]/rhs.m - result.mom[i]/result.m);
    }
    deltal1[0] += r1[1]*p1[2] - r1[2]*p1[1];
    deltal2[0] += r2[1]*p2[2] - r2[2]*p2[1];
    deltal1[1] += r1[2]*p1[0] - r1[0]*p1[2];
    deltal2[1] += r2[2]*p2[0] - r2[0]*p2[2];
    deltal1[2] += r1[0]*p1[1] - r1[1]*p1[0];
    deltal2[2] += r2[0]*p2[1] - r2[1]*p2[0];
    // finally, add the orbital angular momenta
    for (int i=0; i<CH_SPACEDIM; i++) {
      result.angmom[i] += deltal1[i];
      result.angmom[i] += deltal2[i];
    }
  }
#ifdef MERGENEW
  result.newParticle = newParticle && rhs.newParticle;
#endif
  return(result);
}

SinkParticleData&
SinkParticleData::operator+=(const SinkParticleData& rhs) {
  Real x1[CH_SPACEDIM], x2[CH_SPACEDIM], p1[CH_SPACEDIM], p2[CH_SPACEDIM];
  Real m1, m2;
  for (int i=0; i<CH_SPACEDIM; i++) {
    m1 = m;
    m2 = rhs.m;
    x1[i]=pos[i];
    x2[i]=rhs.pos[i];
    p1[i]=mom[i];
    p2[i]=rhs.mom[i];
    pos[i]=pos[i]*m+rhs.pos[i]*rhs.m;
    mom[i]+=rhs.mom[i];
    // put the momentum in the frame of the centroid
    p1[i]=m1*(p1[i]/m1 - mom[i]/(m1+m2));
    p2[i]=m2*(p2[i]/m2 - mom[i]/(m1+m2));
    angmom[i]+=rhs.angmom[i];
  }
  m+=rhs.m;
  for (int i=0; i<CH_SPACEDIM; i++) pos[i] /= m;
  // Add the orbital angular momentum of the two stars about the center
  // of mass to the spin angular momentum of the individual stars
  if (CH_SPACEDIM==3) {
    Real r1[CH_SPACEDIM], r2[CH_SPACEDIM], deltal1[CH_SPACEDIM], deltal2[CH_SPACEDIM];
    for (int i=0; i<CH_SPACEDIM; i++) {
      deltal1[i] = 0.;
      deltal2[i] = 0.;
      r1[i] = x1[i] - pos[i];
      r2[i] = rhs.pos[i] - pos[i];
    }
    deltal1[0] += r1[1]*p1[2] - r1[2]*p1[1];
    deltal2[0] += r2[1]*p2[2] - r2[2]*p2[1];
    deltal1[1] += r1[2]*p1[0] - r1[0]*p1[2];
    deltal2[1] += r2[2]*p2[0] - r2[0]*p2[2];
    deltal1[2] += r1[0]*p1[1] - r1[1]*p1[0];
    deltal2[2] += r2[0]*p2[1] - r2[1]*p2[0];
    // finally, add the orbital angular momenta
    for (int i=0; i<CH_SPACEDIM; i++) {
      angmom[i] += deltal1[i];
      angmom[i] += deltal2[i];
    }
  }
#ifdef MERGENEW
  newParticle = newParticle && rhs.newParticle;
#endif
  return(*this);
}


inline Real* SinkParticleData::momentum() { return(mom); }
inline Real* SinkParticleData::angularMomentum() { return(angmom); }
inline Real SinkParticleData::mass() { return(m); }

void
SinkParticleData::advance(Real dt) {
  //cout << "Entered SinkParticleData::advance" << endl;
  //cout << "SinkParticleData::advance entry: " << pos[0] << endl;
  for (int i=0; i<CH_SPACEDIM; i++) pos[i] += dt*mom[i]/m;
  //cout << "SinkParticleData::advance exit: " << pos[0] << endl;
}
