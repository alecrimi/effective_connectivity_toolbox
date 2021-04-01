window.onload = function() {

  // create and initialize a 3D renderer
  var r = new X.renderer3D();
  r.init();
  r.camera.position = [0, 0, 150];
  
  // create a new X.fibers
  var fibers = new X.fibers();
  // .. associate the TrackVis .TRK file
  fibers.file = 'test.trk';
  fibers.caption = 'The Corpus Callosum:<br>connecting the two hemispheres<br>of the human brain.';
  
  // .. add the fibers
  r.add(fibers);
  
  // showtime!
  r.render();
  
};
