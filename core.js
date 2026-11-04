// ---- Litz physics core (SI units) ----
const MU0 = 4e-7 * Math.PI;
const cadd = (a,b)=>[a[0]+b[0],a[1]+b[1]];
const cmul = (a,b)=>[a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0]];
const cdiv = (a,b)=>{const d=b[0]*b[0]+b[1]*b[1];return [(a[0]*b[0]+a[1]*b[1])/d,(a[1]*b[0]-a[0]*b[1])/d];};
const cabs2 = a=>a[0]*a[0]+a[1]*a[1];
const cconj = a=>[a[0],-a[1]];

// modified Bessel I0, I1 of complex z by power series
function besselI01(z){
  const q = cmul([z[0]/2,z[1]/2],[z[0]/2,z[1]/2]);   // (z/2)^2
  let t0=[1,0], s0=[1,0];
  let t1=[z[0]/2,z[1]/2], s1=[z[0]/2,z[1]/2];
  for(let k=1;k<600;k++){
    t0 = cmul(t0,q); t0=[t0[0]/(k*k),t0[1]/(k*k)];
    t1 = cmul(t1,q); t1=[t1[0]/(k*(k+1)),t1[1]/(k*(k+1))];
    s0=cadd(s0,t0); s1=cadd(s1,t1);
    if(cabs2(t0)<1e-34*cabs2(s0) && cabs2(t1)<1e-34*cabs2(s1) && k>3) break;
  }
  return [s0,s1];
}

function resistivity(T){ return 1.724e-8*(1+0.00393*(T-20)); }   // copper, ohm m

// per-strand AC quantities at frequency f (Hz)
// returns skin factor FR = Rac/Rdc of an isolated strand, and proximity coefficient kp
// kp : loss per metre of strand per (A/m)^2 of RMS transverse field   [W/m / (A/m)^2]
function strandAC(d, f, rho){
  const sigma = 1/rho, w = 2*Math.PI*f, r = d/2;
  const delta = Math.sqrt(2*rho/(w*MU0));
  if(f<=0) return {FR:1,kp:0,delta:Infinity,Gp:1};
  const gamma = [1/delta,1/delta];
  const z = [r/delta, r/delta];
  const [I0,I1] = besselI01(z);
  // skin: Rac/Rdc = Re[ z I0 / (2 I1) ]
  const FR = cdiv(cmul(z,I0),[2*I1[0],2*I1[1]])[0];
  // proximity (uniform transverse field, peak amplitude 1 A/m)
  const gI0 = cmul(gamma,I0);
  const C = cdiv([2*MU0,0],gI0);
  const dI1 = cadd(I0,[ -cdiv(I1,z)[0], -cdiv(I1,z)[1] ]);   // I1'(z)=I0-I1/z
  const prod = cmul(cmul(I1,cconj(gamma)),cconj(dI1));
  const jP = [-prod[1]*cabs2(C), prod[0]*cabs2(C)];            // j*|C|^2*prod
  const Ppk = Math.PI*r*w/(2*MU0)*jP[0];                       // W/m for 1 A/m peak
  const kp = 2*Ppk;                                            // per (A/m rms)^2
  const Plf = sigma*Math.PI*w*w*MU0*MU0*Math.pow(r,4)/8;       // low-freq limit, peak 1 A/m
  return {FR, kp, delta, Gp:Ppk/Plf};
}

// Litz bundle: N strands, strand dia d, bundle dia Db, length l, rms current I
function litz({d,N,Db,len,f,I,T,twist}){
  const rho = resistivity(T), Acu = N*Math.PI*d*d/4;
  const L = len*(twist||1);
  const Rdc = rho*L/Acu;
  const s = strandAC(d,f,rho);
  const Hi2 = I*I/(8*Math.PI*Math.PI*Math.pow(Db/2,2));        // mean-square internal field, A/m rms
  const Pdc = Rdc*I*I;
  const Pskin = Pdc*s.FR;
  const Pprox = N*s.kp*L*Hi2;
  const Ptot = Pskin+Pprox;
  const Rac = Ptot/(I*I);
  // solid round wire of equal copper area
  const dSolid = Math.sqrt(4*Acu/Math.PI);
  const sSolid = strandAC(dSolid,f,rho);
  return {rho,Acu,Rdc,Rac,Pdc,Pskin,Pprox,Ptot,FR:s.FR,Gp:s.Gp,delta:s.delta,
          Hi2,fill:Acu/(Math.PI*Db*Db/4),dSolid,RacSolid:Rdc*sSolid.FR,
          Fr:Rac/Rdc, Rskin:Rdc*s.FR, Rprox:Pprox/(I*I)};
}
if(typeof module!=='undefined') module.exports={litz,strandAC,besselI01,resistivity,MU0};
