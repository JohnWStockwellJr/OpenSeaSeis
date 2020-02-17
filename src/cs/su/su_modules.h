/* Copyright (c) Colorado School of Mines, 2013.*/
/* All rights reserved.                       */

#ifndef SU_MODULES_H
#define SU_MODULES_H

typedef void* (*SUMainMethod) ( void* );

char *sdoc[] = {
"                                                       ",
" GENERIC SU self-documentation                         ",
"                                                       ",
NULL};

const int NUM_SU_MODULES = 234;


std::string SU_MODULE_NAMES[NUM_SU_MODULES] = {
  std::string("sugain"),
  std::string("sudipdivcor"),
  std::string("sucentsamp"),
  std::string("supgc"),
  std::string("sunan"),
  std::string("suweight"),
  std::string("suzero"),
  std::string("sudivcor"),
  std::string("sunormalize"),
  std::string("sucmp"),
  std::string("suquantile"),
  std::string("sumax"),
  std::string("subackus"),
  std::string("sulprime"),
  std::string("sumean"),
  std::string("subackush"),
  std::string("suhistogram"),
  std::string("suattributes"),
  std::string("suacor"),
  std::string("suconv"),
  std::string("surefcon"),
  std::string("suxcor"),
  std::string("suacorfrac"),
  std::string("suunpack1"),
  std::string("supack1"),
  std::string("supack2"),
  std::string("suunpack2"),
  std::string("sustkvel"),
  std::string("segyread"),
  std::string("swapbhed"),
  std::string("segyclean"),
  std::string("suoldtonew"),
  std::string("dt1tosu"),
  std::string("suascii"),
  std::string("suswapbytes"),
  std::string("segyhdrmod"),
  std::string("suintvel"),
  std::string("setbhed"),
  std::string("las2su"),
  std::string("segyhdrs"),
  std::string("segywrite"),
  std::string("sukdmdcs"),
  std::string("sudatumk2ds"),
  std::string("sudatumfd"),
  std::string("sudatumk2dr"),
  std::string("sukdmdcr"),
  std::string("sucddecon"),
  std::string("sufxdecon"),
  std::string("supef"),
  std::string("suphidecon"),
  std::string("sushape"),
  std::string("sudmotx"),
  std::string("sudmofk"),
  std::string("sutihaledmo"),
  std::string("sudmotivz"),
  std::string("sudmofkcw"),
  std::string("sudmovz"),
  std::string("subfilt"),
  std::string("succfilt"),
  std::string("sudipfilt"),
  std::string("sufilter"),
  std::string("sufrac"),
  std::string("sufwatrim"),
  std::string("suk1k2filter"),
  std::string("sukfilter"),
  std::string("sukfrac"),
  std::string("sulfaf"),
  std::string("sumedian"),
  std::string("suphase"),
  std::string("susmgauss2"),
  std::string("sutvband"),
  std::string("su2voxet"),
  std::string("bhedtopar"),
  std::string("su3dchart"),
  std::string("suabshw"),
  std::string("suaddhead"),
  std::string("suahw"),
  std::string("suazimuth"),
  std::string("sucdpbin"),
  std::string("suchart"),
  std::string("suchw"),
  std::string("sucliphead"),
  std::string("sucountkey"),
  std::string("sudumptrace"),
  std::string("sugethw"),
  std::string("suhtmath"),
  std::string("sukeycount"),
  std::string("sulcthw"),
  std::string("sulhead"),
  std::string("supaste"),
  std::string("surandhw"),
  std::string("surange"),
  std::string("susehw"),
  std::string("sushw"),
  std::string("sustrip"),
  std::string("sutrcount"),
  std::string("suutm"),
  std::string("suinterp"),
  std::string("suocext"),
  std::string("suinterpfowler"),
  std::string("sumigtk"),
  std::string("sumigsplit"),
  std::string("sumigffd"),
  std::string("sustolt"),
  std::string("sugazmig"),
  std::string("sumiggbzoan"),
  std::string("sukdmig2d"),
  std::string("sumiggbzo"),
  std::string("sumigpsti"),
  std::string("sumigpreffd"),
  std::string("sumigtopo2d"),
  std::string("sumigpspi"),
  std::string("suinvzco3d"),
  std::string("sumigprefd"),
  std::string("sutifowler"),
  std::string("sumigprepspi"),
  std::string("suktmig2d"),
  std::string("sukdmig3d"),
  std::string("sumigps"),
  std::string("suinvvxzco"),
  std::string("sumigpresp"),
  std::string("sumigfd"),
  std::string("sualford"),
  std::string("sueipofi"),
  std::string("suhrot"),
  std::string("sultt"),
  std::string("supofilt"),
  std::string("supolar"),
  std::string("suaddnoise"),
  std::string("suharlan"),
  std::string("sujitter"),
  std::string("sumix"),
  std::string("suop"),
  std::string("suflip"),
  std::string("suvcat"),
  std::string("suop2"),
  std::string("sufwmix"),
  std::string("supermute"),
  std::string("supickamp"),
  std::string("sufnzero"),
  std::string("sufbpickw"),
  std::string("sucvs4fowler"),
  std::string("sudivstack"),
  std::string("supws"),
  std::string("surecip"),
  std::string("sustack"),
  std::string("sustaticrrs"),
  std::string("suresstat"),
  std::string("sustatic"),
  std::string("surandstat"),
  std::string("suilog"),
  std::string("sureduce"),
  std::string("sutaupnmo"),
  std::string("sushift"),
  std::string("sunmo_a"),
  std::string("suttoz"),
  std::string("sunmo"),
  std::string("sulog"),
  std::string("suresamp"),
  std::string("sutsq"),
  std::string("suztot"),
  std::string("susyncz"),
  std::string("sugassmann"),
  std::string("sukdsyn2d"),
  std::string("suspike"),
  std::string("suimp3d"),
  std::string("suvibro"),
  std::string("surandspike"),
  std::string("suimp2d"),
  std::string("susynlvfti"),
  std::string("suwaveform"),
  std::string("sunhmospike"),
  std::string("susynlv"),
  std::string("sufctanismod"),
  std::string("sugoupillaud"),
  std::string("susynlvcw"),
  std::string("sufdmod2_pml"),
  std::string("susynvxzcs"),
  std::string("sufdmod2"),
  std::string("susynvxz"),
  std::string("suaddevent"),
  std::string("suimpedance"),
  std::string("sunull"),
  std::string("suplane"),
  std::string("sugoupillaudpo"),
  std::string("sufdmod1"),
  std::string("suea2df"),
  std::string("suramp"),
  std::string("sutxtaper"),
  std::string("sutaper"),
  std::string("sugausstaper"),
  std::string("suamp"),
  std::string("suanalytic"),
  std::string("succepstrum"),
  std::string("succwt"),
  std::string("sucepstrum"),
  std::string("suclogfft"),
  std::string("sucwt"),
  std::string("sufft"),
  std::string("sugabor"),
  std::string("suhilb"),
  std::string("suiclogfft"),
  std::string("suifft"),
  std::string("suphasevel"),
  std::string("suradon"),
  std::string("suslowft"),
  std::string("suslowift"),
  std::string("suspecfk"),
  std::string("suspecfx"),
  std::string("suspeck1k2"),
  std::string("sutaup"),
  std::string("suwfft"),
  std::string("sutivel"),
  std::string("suvel2df"),
  std::string("surelan"),
  std::string("suvelan_nccs"),
  std::string("surelanan"),
  std::string("suvelan_nsel"),
  std::string("suvelan_uccs"),
  std::string("suvelan"),
  std::string("suvelan_usel"),
  std::string("suwellrf"),
  std::string("sukill"),
  std::string("sugprfb"),
  std::string("supad"),
  std::string("sucommand"),
  std::string("susorty"),
  std::string("suwindpoly"),
  std::string("susort"),
  std::string("suputgthr"),
  std::string("susplit"),
  std::string("sumute"),
  std::string("sumixgathers"),
  std::string("suwind")
};
namespace sugain {
  void* main_sugain( void* args );
}
namespace sudipdivcor {
  void* main_sudipdivcor( void* args );
}
namespace sucentsamp {
  void* main_sucentsamp( void* args );
}
namespace supgc {
  void* main_supgc( void* args );
}
namespace sunan {
  void* main_sunan( void* args );
}
namespace suweight {
  void* main_suweight( void* args );
}
namespace suzero {
  void* main_suzero( void* args );
}
namespace sudivcor {
  void* main_sudivcor( void* args );
}
namespace sunormalize {
  void* main_sunormalize( void* args );
}
namespace sucmp {
  void* main_sucmp( void* args );
}
namespace suquantile {
  void* main_suquantile( void* args );
}
namespace sumax {
  void* main_sumax( void* args );
}
namespace subackus {
  void* main_subackus( void* args );
}
namespace sulprime {
  void* main_sulprime( void* args );
}
namespace sumean {
  void* main_sumean( void* args );
}
namespace subackush {
  void* main_subackush( void* args );
}
namespace suhistogram {
  void* main_suhistogram( void* args );
}
namespace suattributes {
  void* main_suattributes( void* args );
}
namespace suacor {
  void* main_suacor( void* args );
}
namespace suconv {
  void* main_suconv( void* args );
}
namespace surefcon {
  void* main_surefcon( void* args );
}
namespace suxcor {
  void* main_suxcor( void* args );
}
namespace suacorfrac {
  void* main_suacorfrac( void* args );
}
namespace suunpack1 {
  void* main_suunpack1( void* args );
}
namespace supack1 {
  void* main_supack1( void* args );
}
namespace supack2 {
  void* main_supack2( void* args );
}
namespace suunpack2 {
  void* main_suunpack2( void* args );
}
namespace sustkvel {
  void* main_sustkvel( void* args );
}
namespace segyread {
  void* main_segyread( void* args );
}
namespace swapbhed {
  void* main_swapbhed( void* args );
}
namespace segyclean {
  void* main_segyclean( void* args );
}
namespace suoldtonew {
  void* main_suoldtonew( void* args );
}
namespace dt1tosu {
  void* main_dt1tosu( void* args );
}
namespace suascii {
  void* main_suascii( void* args );
}
namespace suswapbytes {
  void* main_suswapbytes( void* args );
}
namespace segyhdrmod {
  void* main_segyhdrmod( void* args );
}
namespace suintvel {
  void* main_suintvel( void* args );
}
namespace setbhed {
  void* main_setbhed( void* args );
}
namespace las2su {
  void* main_las2su( void* args );
}
namespace segyhdrs {
  void* main_segyhdrs( void* args );
}
namespace segywrite {
  void* main_segywrite( void* args );
}
namespace sukdmdcs {
  void* main_sukdmdcs( void* args );
}
namespace sudatumk2ds {
  void* main_sudatumk2ds( void* args );
}
namespace sudatumfd {
  void* main_sudatumfd( void* args );
}
namespace sudatumk2dr {
  void* main_sudatumk2dr( void* args );
}
namespace sukdmdcr {
  void* main_sukdmdcr( void* args );
}
namespace sucddecon {
  void* main_sucddecon( void* args );
}
namespace sufxdecon {
  void* main_sufxdecon( void* args );
}
namespace supef {
  void* main_supef( void* args );
}
namespace suphidecon {
  void* main_suphidecon( void* args );
}
namespace sushape {
  void* main_sushape( void* args );
}
namespace sudmotx {
  void* main_sudmotx( void* args );
}
namespace sudmofk {
  void* main_sudmofk( void* args );
}
namespace sutihaledmo {
  void* main_sutihaledmo( void* args );
}
namespace sudmotivz {
  void* main_sudmotivz( void* args );
}
namespace sudmofkcw {
  void* main_sudmofkcw( void* args );
}
namespace sudmovz {
  void* main_sudmovz( void* args );
}
namespace subfilt {
  void* main_subfilt( void* args );
}
namespace succfilt {
  void* main_succfilt( void* args );
}
namespace sudipfilt {
  void* main_sudipfilt( void* args );
}
namespace sufilter {
  void* main_sufilter( void* args );
}
namespace sufrac {
  void* main_sufrac( void* args );
}
namespace sufwatrim {
  void* main_sufwatrim( void* args );
}
namespace suk1k2filter {
  void* main_suk1k2filter( void* args );
}
namespace sukfilter {
  void* main_sukfilter( void* args );
}
namespace sukfrac {
  void* main_sukfrac( void* args );
}
namespace sulfaf {
  void* main_sulfaf( void* args );
}
namespace sumedian {
  void* main_sumedian( void* args );
}
namespace suphase {
  void* main_suphase( void* args );
}
namespace susmgauss2 {
  void* main_susmgauss2( void* args );
}
namespace sutvband {
  void* main_sutvband( void* args );
}
namespace su2voxet {
  void* main_su2voxet( void* args );
}
namespace bhedtopar {
  void* main_bhedtopar( void* args );
}
namespace su3dchart {
  void* main_su3dchart( void* args );
}
namespace suabshw {
  void* main_suabshw( void* args );
}
namespace suaddhead {
  void* main_suaddhead( void* args );
}
namespace suahw {
  void* main_suahw( void* args );
}
namespace suazimuth {
  void* main_suazimuth( void* args );
}
namespace sucdpbin {
  void* main_sucdpbin( void* args );
}
namespace suchart {
  void* main_suchart( void* args );
}
namespace suchw {
  void* main_suchw( void* args );
}
namespace sucliphead {
  void* main_sucliphead( void* args );
}
namespace sucountkey {
  void* main_sucountkey( void* args );
}
namespace sudumptrace {
  void* main_sudumptrace( void* args );
}
namespace sugethw {
  void* main_sugethw( void* args );
}
namespace suhtmath {
  void* main_suhtmath( void* args );
}
namespace sukeycount {
  void* main_sukeycount( void* args );
}
namespace sulcthw {
  void* main_sulcthw( void* args );
}
namespace sulhead {
  void* main_sulhead( void* args );
}
namespace supaste {
  void* main_supaste( void* args );
}
namespace surandhw {
  void* main_surandhw( void* args );
}
namespace surange {
  void* main_surange( void* args );
}
namespace susehw {
  void* main_susehw( void* args );
}
namespace sushw {
  void* main_sushw( void* args );
}
namespace sustrip {
  void* main_sustrip( void* args );
}
namespace sutrcount {
  void* main_sutrcount( void* args );
}
namespace suutm {
  void* main_suutm( void* args );
}
namespace suinterp {
  void* main_suinterp( void* args );
}
namespace suocext {
  void* main_suocext( void* args );
}
namespace suinterpfowler {
  void* main_suinterpfowler( void* args );
}
namespace sumigtk {
  void* main_sumigtk( void* args );
}
namespace sumigsplit {
  void* main_sumigsplit( void* args );
}
namespace sumigffd {
  void* main_sumigffd( void* args );
}
namespace sustolt {
  void* main_sustolt( void* args );
}
namespace sugazmig {
  void* main_sugazmig( void* args );
}
namespace sumiggbzoan {
  void* main_sumiggbzoan( void* args );
}
namespace sukdmig2d {
  void* main_sukdmig2d( void* args );
}
namespace sumiggbzo {
  void* main_sumiggbzo( void* args );
}
namespace sumigpsti {
  void* main_sumigpsti( void* args );
}
namespace sumigpreffd {
  void* main_sumigpreffd( void* args );
}
namespace sumigtopo2d {
  void* main_sumigtopo2d( void* args );
}
namespace sumigpspi {
  void* main_sumigpspi( void* args );
}
namespace suinvzco3d {
  void* main_suinvzco3d( void* args );
}
namespace sumigprefd {
  void* main_sumigprefd( void* args );
}
namespace sutifowler {
  void* main_sutifowler( void* args );
}
namespace sumigprepspi {
  void* main_sumigprepspi( void* args );
}
namespace suktmig2d {
  void* main_suktmig2d( void* args );
}
namespace sukdmig3d {
  void* main_sukdmig3d( void* args );
}
namespace sumigps {
  void* main_sumigps( void* args );
}
namespace suinvvxzco {
  void* main_suinvvxzco( void* args );
}
namespace sumigpresp {
  void* main_sumigpresp( void* args );
}
namespace sumigfd {
  void* main_sumigfd( void* args );
}
namespace sualford {
  void* main_sualford( void* args );
}
namespace sueipofi {
  void* main_sueipofi( void* args );
}
namespace suhrot {
  void* main_suhrot( void* args );
}
namespace sultt {
  void* main_sultt( void* args );
}
namespace supofilt {
  void* main_supofilt( void* args );
}
namespace supolar {
  void* main_supolar( void* args );
}
namespace suaddnoise {
  void* main_suaddnoise( void* args );
}
namespace suharlan {
  void* main_suharlan( void* args );
}
namespace sujitter {
  void* main_sujitter( void* args );
}
namespace sumix {
  void* main_sumix( void* args );
}
namespace suop {
  void* main_suop( void* args );
}
namespace suflip {
  void* main_suflip( void* args );
}
namespace suvcat {
  void* main_suvcat( void* args );
}
namespace suop2 {
  void* main_suop2( void* args );
}
namespace sufwmix {
  void* main_sufwmix( void* args );
}
namespace supermute {
  void* main_supermute( void* args );
}
namespace supickamp {
  void* main_supickamp( void* args );
}
namespace sufnzero {
  void* main_sufnzero( void* args );
}
namespace sufbpickw {
  void* main_sufbpickw( void* args );
}
namespace sucvs4fowler {
  void* main_sucvs4fowler( void* args );
}
namespace sudivstack {
  void* main_sudivstack( void* args );
}
namespace supws {
  void* main_supws( void* args );
}
namespace surecip {
  void* main_surecip( void* args );
}
namespace sustack {
  void* main_sustack( void* args );
}
namespace sustaticrrs {
  void* main_sustaticrrs( void* args );
}
namespace suresstat {
  void* main_suresstat( void* args );
}
namespace sustatic {
  void* main_sustatic( void* args );
}
namespace surandstat {
  void* main_surandstat( void* args );
}
namespace suilog {
  void* main_suilog( void* args );
}
namespace sureduce {
  void* main_sureduce( void* args );
}
namespace sutaupnmo {
  void* main_sutaupnmo( void* args );
}
namespace sushift {
  void* main_sushift( void* args );
}
namespace sunmo_a {
  void* main_sunmo_a( void* args );
}
namespace suttoz {
  void* main_suttoz( void* args );
}
namespace sunmo {
  void* main_sunmo( void* args );
}
namespace sulog {
  void* main_sulog( void* args );
}
namespace suresamp {
  void* main_suresamp( void* args );
}
namespace sutsq {
  void* main_sutsq( void* args );
}
namespace suztot {
  void* main_suztot( void* args );
}
namespace susyncz {
  void* main_susyncz( void* args );
}
namespace sugassmann {
  void* main_sugassmann( void* args );
}
namespace sukdsyn2d {
  void* main_sukdsyn2d( void* args );
}
namespace suspike {
  void* main_suspike( void* args );
}
namespace suimp3d {
  void* main_suimp3d( void* args );
}
namespace suvibro {
  void* main_suvibro( void* args );
}
namespace surandspike {
  void* main_surandspike( void* args );
}
namespace suimp2d {
  void* main_suimp2d( void* args );
}
namespace susynlvfti {
  void* main_susynlvfti( void* args );
}
namespace suwaveform {
  void* main_suwaveform( void* args );
}
namespace sunhmospike {
  void* main_sunhmospike( void* args );
}
namespace susynlv {
  void* main_susynlv( void* args );
}
namespace sufctanismod {
  void* main_sufctanismod( void* args );
}
namespace sugoupillaud {
  void* main_sugoupillaud( void* args );
}
namespace susynlvcw {
  void* main_susynlvcw( void* args );
}
namespace sufdmod2_pml {
  void* main_sufdmod2_pml( void* args );
}
namespace susynvxzcs {
  void* main_susynvxzcs( void* args );
}
namespace sufdmod2 {
  void* main_sufdmod2( void* args );
}
namespace susynvxz {
  void* main_susynvxz( void* args );
}
namespace suaddevent {
  void* main_suaddevent( void* args );
}
namespace suimpedance {
  void* main_suimpedance( void* args );
}
namespace sunull {
  void* main_sunull( void* args );
}
namespace suplane {
  void* main_suplane( void* args );
}
namespace sugoupillaudpo {
  void* main_sugoupillaudpo( void* args );
}
namespace sufdmod1 {
  void* main_sufdmod1( void* args );
}
namespace suea2df {
  void* main_suea2df( void* args );
}
namespace suramp {
  void* main_suramp( void* args );
}
namespace sutxtaper {
  void* main_sutxtaper( void* args );
}
namespace sutaper {
  void* main_sutaper( void* args );
}
namespace sugausstaper {
  void* main_sugausstaper( void* args );
}
namespace suamp {
  void* main_suamp( void* args );
}
namespace suanalytic {
  void* main_suanalytic( void* args );
}
namespace succepstrum {
  void* main_succepstrum( void* args );
}
namespace succwt {
  void* main_succwt( void* args );
}
namespace sucepstrum {
  void* main_sucepstrum( void* args );
}
namespace suclogfft {
  void* main_suclogfft( void* args );
}
namespace sucwt {
  void* main_sucwt( void* args );
}
namespace sufft {
  void* main_sufft( void* args );
}
namespace sugabor {
  void* main_sugabor( void* args );
}
namespace suhilb {
  void* main_suhilb( void* args );
}
namespace suiclogfft {
  void* main_suiclogfft( void* args );
}
namespace suifft {
  void* main_suifft( void* args );
}
namespace suphasevel {
  void* main_suphasevel( void* args );
}
namespace suradon {
  void* main_suradon( void* args );
}
namespace suslowft {
  void* main_suslowft( void* args );
}
namespace suslowift {
  void* main_suslowift( void* args );
}
namespace suspecfk {
  void* main_suspecfk( void* args );
}
namespace suspecfx {
  void* main_suspecfx( void* args );
}
namespace suspeck1k2 {
  void* main_suspeck1k2( void* args );
}
namespace sutaup {
  void* main_sutaup( void* args );
}
namespace suwfft {
  void* main_suwfft( void* args );
}
namespace sutivel {
  void* main_sutivel( void* args );
}
namespace suvel2df {
  void* main_suvel2df( void* args );
}
namespace surelan {
  void* main_surelan( void* args );
}
namespace suvelan_nccs {
  void* main_suvelan_nccs( void* args );
}
namespace surelanan {
  void* main_surelanan( void* args );
}
namespace suvelan_nsel {
  void* main_suvelan_nsel( void* args );
}
namespace suvelan_uccs {
  void* main_suvelan_uccs( void* args );
}
namespace suvelan {
  void* main_suvelan( void* args );
}
namespace suvelan_usel {
  void* main_suvelan_usel( void* args );
}
namespace suwellrf {
  void* main_suwellrf( void* args );
}
namespace sukill {
  void* main_sukill( void* args );
}
namespace sugprfb {
  void* main_sugprfb( void* args );
}
namespace supad {
  void* main_supad( void* args );
}
namespace sucommand {
  void* main_sucommand( void* args );
}
namespace susorty {
  void* main_susorty( void* args );
}
namespace suwindpoly {
  void* main_suwindpoly( void* args );
}
namespace susort {
  void* main_susort( void* args );
}
namespace suputgthr {
  void* main_suputgthr( void* args );
}
namespace susplit {
  void* main_susplit( void* args );
}
namespace sumute {
  void* main_sumute( void* args );
}
namespace sumixgathers {
  void* main_sumixgathers( void* args );
}
namespace suwind {
  void* main_suwind( void* args );
}

void*(*SU_MAIN_METHODS[NUM_SU_MODULES])( void* args ) = {
  sugain::main_sugain,
  sudipdivcor::main_sudipdivcor,
  sucentsamp::main_sucentsamp,
  supgc::main_supgc,
  sunan::main_sunan,
  suweight::main_suweight,
  suzero::main_suzero,
  sudivcor::main_sudivcor,
  sunormalize::main_sunormalize,
  sucmp::main_sucmp,
  suquantile::main_suquantile,
  sumax::main_sumax,
  subackus::main_subackus,
  sulprime::main_sulprime,
  sumean::main_sumean,
  subackush::main_subackush,
  suhistogram::main_suhistogram,
  suattributes::main_suattributes,
  suacor::main_suacor,
  suconv::main_suconv,
  surefcon::main_surefcon,
  suxcor::main_suxcor,
  suacorfrac::main_suacorfrac,
  suunpack1::main_suunpack1,
  supack1::main_supack1,
  supack2::main_supack2,
  suunpack2::main_suunpack2,
  sustkvel::main_sustkvel,
  segyread::main_segyread,
  swapbhed::main_swapbhed,
  segyclean::main_segyclean,
  suoldtonew::main_suoldtonew,
  dt1tosu::main_dt1tosu,
  suascii::main_suascii,
  suswapbytes::main_suswapbytes,
  segyhdrmod::main_segyhdrmod,
  suintvel::main_suintvel,
  setbhed::main_setbhed,
  las2su::main_las2su,
  segyhdrs::main_segyhdrs,
  segywrite::main_segywrite,
  sukdmdcs::main_sukdmdcs,
  sudatumk2ds::main_sudatumk2ds,
  sudatumfd::main_sudatumfd,
  sudatumk2dr::main_sudatumk2dr,
  sukdmdcr::main_sukdmdcr,
  sucddecon::main_sucddecon,
  sufxdecon::main_sufxdecon,
  supef::main_supef,
  suphidecon::main_suphidecon,
  sushape::main_sushape,
  sudmotx::main_sudmotx,
  sudmofk::main_sudmofk,
  sutihaledmo::main_sutihaledmo,
  sudmotivz::main_sudmotivz,
  sudmofkcw::main_sudmofkcw,
  sudmovz::main_sudmovz,
  subfilt::main_subfilt,
  succfilt::main_succfilt,
  sudipfilt::main_sudipfilt,
  sufilter::main_sufilter,
  sufrac::main_sufrac,
  sufwatrim::main_sufwatrim,
  suk1k2filter::main_suk1k2filter,
  sukfilter::main_sukfilter,
  sukfrac::main_sukfrac,
  sulfaf::main_sulfaf,
  sumedian::main_sumedian,
  suphase::main_suphase,
  susmgauss2::main_susmgauss2,
  sutvband::main_sutvband,
  su2voxet::main_su2voxet,
  bhedtopar::main_bhedtopar,
  su3dchart::main_su3dchart,
  suabshw::main_suabshw,
  suaddhead::main_suaddhead,
  suahw::main_suahw,
  suazimuth::main_suazimuth,
  sucdpbin::main_sucdpbin,
  suchart::main_suchart,
  suchw::main_suchw,
  sucliphead::main_sucliphead,
  sucountkey::main_sucountkey,
  sudumptrace::main_sudumptrace,
  sugethw::main_sugethw,
  suhtmath::main_suhtmath,
  sukeycount::main_sukeycount,
  sulcthw::main_sulcthw,
  sulhead::main_sulhead,
  supaste::main_supaste,
  surandhw::main_surandhw,
  surange::main_surange,
  susehw::main_susehw,
  sushw::main_sushw,
  sustrip::main_sustrip,
  sutrcount::main_sutrcount,
  suutm::main_suutm,
  suinterp::main_suinterp,
  suocext::main_suocext,
  suinterpfowler::main_suinterpfowler,
  sumigtk::main_sumigtk,
  sumigsplit::main_sumigsplit,
  sumigffd::main_sumigffd,
  sustolt::main_sustolt,
  sugazmig::main_sugazmig,
  sumiggbzoan::main_sumiggbzoan,
  sukdmig2d::main_sukdmig2d,
  sumiggbzo::main_sumiggbzo,
  sumigpsti::main_sumigpsti,
  sumigpreffd::main_sumigpreffd,
  sumigtopo2d::main_sumigtopo2d,
  sumigpspi::main_sumigpspi,
  suinvzco3d::main_suinvzco3d,
  sumigprefd::main_sumigprefd,
  sutifowler::main_sutifowler,
  sumigprepspi::main_sumigprepspi,
  suktmig2d::main_suktmig2d,
  sukdmig3d::main_sukdmig3d,
  sumigps::main_sumigps,
  suinvvxzco::main_suinvvxzco,
  sumigpresp::main_sumigpresp,
  sumigfd::main_sumigfd,
  sualford::main_sualford,
  sueipofi::main_sueipofi,
  suhrot::main_suhrot,
  sultt::main_sultt,
  supofilt::main_supofilt,
  supolar::main_supolar,
  suaddnoise::main_suaddnoise,
  suharlan::main_suharlan,
  sujitter::main_sujitter,
  sumix::main_sumix,
  suop::main_suop,
  suflip::main_suflip,
  suvcat::main_suvcat,
  suop2::main_suop2,
  sufwmix::main_sufwmix,
  supermute::main_supermute,
  supickamp::main_supickamp,
  sufnzero::main_sufnzero,
  sufbpickw::main_sufbpickw,
  sucvs4fowler::main_sucvs4fowler,
  sudivstack::main_sudivstack,
  supws::main_supws,
  surecip::main_surecip,
  sustack::main_sustack,
  sustaticrrs::main_sustaticrrs,
  suresstat::main_suresstat,
  sustatic::main_sustatic,
  surandstat::main_surandstat,
  suilog::main_suilog,
  sureduce::main_sureduce,
  sutaupnmo::main_sutaupnmo,
  sushift::main_sushift,
  sunmo_a::main_sunmo_a,
  suttoz::main_suttoz,
  sunmo::main_sunmo,
  sulog::main_sulog,
  suresamp::main_suresamp,
  sutsq::main_sutsq,
  suztot::main_suztot,
  susyncz::main_susyncz,
  sugassmann::main_sugassmann,
  sukdsyn2d::main_sukdsyn2d,
  suspike::main_suspike,
  suimp3d::main_suimp3d,
  suvibro::main_suvibro,
  surandspike::main_surandspike,
  suimp2d::main_suimp2d,
  susynlvfti::main_susynlvfti,
  suwaveform::main_suwaveform,
  sunhmospike::main_sunhmospike,
  susynlv::main_susynlv,
  sufctanismod::main_sufctanismod,
  sugoupillaud::main_sugoupillaud,
  susynlvcw::main_susynlvcw,
  sufdmod2_pml::main_sufdmod2_pml,
  susynvxzcs::main_susynvxzcs,
  sufdmod2::main_sufdmod2,
  susynvxz::main_susynvxz,
  suaddevent::main_suaddevent,
  suimpedance::main_suimpedance,
  sunull::main_sunull,
  suplane::main_suplane,
  sugoupillaudpo::main_sugoupillaudpo,
  sufdmod1::main_sufdmod1,
  suea2df::main_suea2df,
  suramp::main_suramp,
  sutxtaper::main_sutxtaper,
  sutaper::main_sutaper,
  sugausstaper::main_sugausstaper,
  suamp::main_suamp,
  suanalytic::main_suanalytic,
  succepstrum::main_succepstrum,
  succwt::main_succwt,
  sucepstrum::main_sucepstrum,
  suclogfft::main_suclogfft,
  sucwt::main_sucwt,
  sufft::main_sufft,
  sugabor::main_sugabor,
  suhilb::main_suhilb,
  suiclogfft::main_suiclogfft,
  suifft::main_suifft,
  suphasevel::main_suphasevel,
  suradon::main_suradon,
  suslowft::main_suslowft,
  suslowift::main_suslowift,
  suspecfk::main_suspecfk,
  suspecfx::main_suspecfx,
  suspeck1k2::main_suspeck1k2,
  sutaup::main_sutaup,
  suwfft::main_suwfft,
  sutivel::main_sutivel,
  suvel2df::main_suvel2df,
  surelan::main_surelan,
  suvelan_nccs::main_suvelan_nccs,
  surelanan::main_surelanan,
  suvelan_nsel::main_suvelan_nsel,
  suvelan_uccs::main_suvelan_uccs,
  suvelan::main_suvelan,
  suvelan_usel::main_suvelan_usel,
  suwellrf::main_suwellrf,
  sukill::main_sukill,
  sugprfb::main_sugprfb,
  supad::main_supad,
  sucommand::main_sucommand,
  susorty::main_susorty,
  suwindpoly::main_suwindpoly,
  susort::main_susort,
  suputgthr::main_suputgthr,
  susplit::main_susplit,
  sumute::main_sumute,
  sumixgathers::main_sumixgathers,
  suwind::main_suwind
};
#endif
