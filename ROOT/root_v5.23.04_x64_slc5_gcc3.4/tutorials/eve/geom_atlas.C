// @(#)root/eve:$Id: geom_atlas.C 27557 2009-02-20 18:22:42Z matevz $
// Author: Matevz Tadel

// Shows ATLAS geometry.

void geom_atlas()
{
   TEveManager::Create();

   gGeoManager = gEve->GetGeometry("http://root.cern.ch/files/atlas.root");

   TGeoNode* node1 = gGeoManager->GetTopVolume()->FindNode("INNE_1");
   TEveGeoTopNode* inn = new TEveGeoTopNode(gGeoManager, node1);
   gEve->AddGlobalElement(inn);

   TGeoNode* node2 = gGeoManager->GetTopVolume()->FindNode("CENT_1");
   TEveGeoTopNode* cnt = new TEveGeoTopNode(gGeoManager, node2);
   gEve->AddGlobalElement(cnt);

   TGeoNode* node3 = gGeoManager->GetTopVolume()->FindNode("OUTE_1");
   TEveGeoTopNode* out = new TEveGeoTopNode(gGeoManager, node3);
   gEve->AddGlobalElement(out);

   gEve->FullRedraw3D(kTRUE);

   // EClipType not exported to CINT (see TGLUtil.h):
   // 0 - no clip, 1 - clip plane, 2 - clip box
   TGLViewer *v = gEve->GetDefaultGLViewer();
   v->GetClipSet()->SetClipType(1);
   v->RefreshPadEditor(v);

   v->CurrentCamera().RotateRad(-.7, 0.5);
   v->DoDraw();
}
