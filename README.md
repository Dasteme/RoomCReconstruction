# RoomCReconstruction

RoomCReconstruction is a project which is able to reconstruct the planimetry of non-empty rooms.<br>
For this, it takes a point cloud as input and outputs a mesh.<br>
It is divided into three steps:<br>

 1. A Clustering algorithm which is able to separate the point cloud into "clusters" (points that lie in the same plane)
 2. An algorithm called "triangle-finding" which finds corners in the room.
 3. An algorithm called "triangle-linking" which connects the triangles (corners) to a closed room.

The linked triangles are then converted to a mesh.
