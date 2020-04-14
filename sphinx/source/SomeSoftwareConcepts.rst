==================================
Some Solar Corona Raytrace Concept
==================================


Calling sequence of a typical raytracing
========================================

* routine name (filename without extension) : description
* rtthread (rtthread) : main wrapper
* Scene (scene) : sets the scene: camera FOV, density model, attitude, physics, ...
* computeImagebyChunk (scene) : multi-threaded ray-tracing engine
* losintegchunk (scene) : compute a sub frame of the full FOV
* losinteg (scene) : ray-trace one pixel
* computeRadiation (physicsbase (virtual: physics dependent)) : computes the scattered radiation level based on the requested physics
* Density (CModelBase (virtual: model dependent)) : computes the electron or dust density, based on the requested model
