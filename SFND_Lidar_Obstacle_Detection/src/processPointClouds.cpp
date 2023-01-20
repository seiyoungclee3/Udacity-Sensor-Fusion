// PCL lib Functions for processing point clouds 

#include "processPointClouds.h"


//constructor:
template<typename PointT>
ProcessPointClouds<PointT>::ProcessPointClouds() {}


//de-constructor:
template<typename PointT>
ProcessPointClouds<PointT>::~ProcessPointClouds() {}


template<typename PointT>
void ProcessPointClouds<PointT>::numPoints(typename pcl::PointCloud<PointT>::Ptr cloud)
{
    std::cout << cloud->points.size() << std::endl;
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::FilterCloud(typename pcl::PointCloud<PointT>::Ptr cloud, float filterRes, Eigen::Vector4f minPoint, Eigen::Vector4f maxPoint)
{

    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();

    // TODO:: Fill in the function to do voxel grid point reduction and region based filtering
    pcl::VoxelGrid<PointT> vg; 
    typename pcl::PointCloud<PointT>::Ptr cloudFiltered (new pcl::PointCloud<PointT>); 
    vg.setInputCloud(cloud); 
    vg.setLeafSize(filterRes,filterRes,filterRes); 
    vg.filter(*cloudFiltered); 

    typename pcl::PointCloud<PointT>::Ptr cloudRegion (new pcl::PointCloud<PointT>); 

    pcl::CropBox<PointT> region(true); 
    region.setMin(minPoint); 
    region.setMax(maxPoint); 
    region.setInputCloud(cloudFiltered); 
    region.filter(*cloudRegion); 

    std::vector<int> indices; 

    pcl::CropBox<PointT> roof(true);
    roof.setMin(Eigen::Vector4f (-1.5,-1.7,-1,1)); 
    roof.setMax(Eigen::Vector4f (2.6,1.7,-.4,1));
    roof.setInputCloud(cloudRegion); 
    roof.filter(indices); 
    
    pcl::PointIndices::Ptr inliers {new pcl::PointIndices}; 
    for(int point : indices)
        inliers->indices.push_back(point); 
    pcl::ExtractIndices<PointT> extract; 
    extract.setInputCloud (cloudRegion); 
    extract.setIndices (inliers); 
    extract.setNegative (true); 
    extract.filter (*cloudRegion); 

    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "filtering took " << elapsedTime.count() << " milliseconds" << std::endl;

    return cloud;

}

//Extending RANSAC to Planes 
template <typename PointT>
std::unordered_set<int> ProcessPointClouds<PointT>::Ransac3DPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceTol)
{
	std::unordered_set<int> inliersResult;
	srand(time(NULL));
	// For max iterations 
    while(maxIterations--)
	{
		// Randomly sample subset 
		std::unordered_set<int> inliers; 
		while (inliers.size()<3)
			inliers.insert(rand()%(cloud->points.size())); 
		float x1,y1,z1,x2,y2,z2,x3,y3,z3; 

        auto itr = inliers.begin(); 
		x1 = cloud->points[*itr].x; 
		y1 = cloud->points[*itr].y; 
		z1 = cloud->points[*itr].z; 
		itr++; 
		x2 = cloud->points[*itr].x; 
		y2 = cloud->points[*itr].y; 
		z2 = cloud->points[*itr].z; 
		itr++; 
		x3 = cloud->points[*itr].x; 
		y3 = cloud->points[*itr].y; 
		z3 = cloud->points[*itr].z; 

		float A = (y2-y1)*(z3-z1)-(z2-z1)*(y3-y1);
		float B = (z2-z1)*(x3-x1)-(x2-x1)*(z3-z1);
		float C = (x2-x1)*(y3-y1)-(y2-y1)*(x3-x1);
		float D = -(A*x1 +B*y1+C*z1); 
		for(int index =0; index <cloud->points.size();index++){
			if(inliers.count(index)>0)
				continue; 
			
			float x4 = cloud->points[index].x; 
			float y4 = cloud->points[index].y; 
			float z4 = cloud->points[index].z; 

			float d = fabs(A*x4+B*y4+C*z4+D)/sqrt(A*A+B*B+C*C); 
			if(d <= distanceTol)
			{
				inliers.insert(index); 
			}
				
		}

        if(inliers.size()>inliersResult.size()){
			inliersResult=inliers; 
		}
	}
	
	return inliersResult; 

}

template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::SeparateClouds(pcl::PointIndices::Ptr inliers, typename pcl::PointCloud<PointT>::Ptr cloud) 
{
  // TODO: Create two new point clouds, one cloud with obstacles and other with segmented plane
    typename pcl::PointCloud<PointT>::Ptr obstCloud (new pcl::PointCloud<PointT> ());  //obstacle (non-plane cloud )
    typename pcl::PointCloud<PointT>::Ptr planeCloud (new pcl::PointCloud<PointT> ()); //road points (plane cloud)
    
    for(int index: inliers->indices)
        planeCloud->points.push_back(cloud->points[index]); 
    // to generate the obstacle cloud, one way to use PCL to do this is to use an extract object, which subtracts the plane cloud from the input cloud. 
    pcl::ExtractIndices<PointT> extract; 
    extract.setInputCloud (cloud);
    extract.setIndices (inliers);
    extract.setNegative (true);
    extract.filter (*obstCloud); 
    
    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> segResult(obstCloud, planeCloud);
    return segResult;
}


template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::SegmentPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceThreshold)
{
    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();
	
    // lecture: segmenting the plane with PCL 
    // instruction: implement the segmentPlane that will return a  std::pair holding point
    //              cloud types  
    pcl::SACSegmentation<PointT> seg;
    pcl::ModelCoefficients::Ptr coefficients {new pcl::ModelCoefficients ()}; 
    pcl::PointIndices::Ptr inliers {new pcl::PointIndices}; //contains indices of all points 
   
    seg.setOptimizeCoefficients(true); 
    seg.setModelType(pcl::SACMODEL_PLANE); 
    seg.setMethodType(pcl::SAC_RANSAC);
    seg.setMaxIterations(maxIterations);
    seg.setDistanceThreshold(distanceThreshold); 

    //segment the largest planar component from the input cloud 
    seg.setInputCloud(cloud); 
    seg.segment(*inliers, *coefficients);
    if (inliers->indices.size()==0)
    {
        std::cout <<"could not estimate a planar model for the given datasets." <<std::endl;
    }
     
    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> segResult = SeparateClouds(inliers,cloud);
 
    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "plane segmentation took " << elapsedTime.count() << " milliseconds" << std::endl;
      
    return segResult;
}

template<typename PointT>
std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::Segment3DPlane(typename pcl::PointCloud<PointT>::Ptr cloud, int maxIterations, float distanceThreshold)
{
    // Time segmentation process
    auto startTime = std::chrono::steady_clock::now();
    
    pcl::PointIndices::Ptr inliers {new pcl::PointIndices}; //contains indices of all points 

    std::unordered_set<int> inliersRes = Ransac3DPlane(cloud, maxIterations, distanceThreshold);
    for (int i:inliersRes){
        inliers->indices.push_back(i); 
    }

    if (inliers->indices.size()==0)
    {
        std::cout <<"could not estimate a planar model for the given datasets." <<std::endl;
    }
     
    std::pair<typename pcl::PointCloud<PointT>::Ptr, typename pcl::PointCloud<PointT>::Ptr> segResult = SeparateClouds(inliers,cloud);
 
    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "plane segmentation took " << elapsedTime.count() << " milliseconds" << std::endl;
      
    return segResult;
}

template<typename PointT>
std::vector<typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::Clustering(typename pcl::PointCloud<PointT>::Ptr cloud, float clusterTolerance, int minSize, int maxSize)
{
    // Time clustering process
    auto startTime = std::chrono::steady_clock::now();

    std::vector<typename pcl::PointCloud<PointT>::Ptr> clusters;

    // TODO:: Fill in the function to perform euclidean clustering to group detected obstacles
    // creating the KDTree object for the search method of the extraction 
    typename pcl::search::KdTree<PointT>::Ptr tree(new pcl::search::KdTree<PointT>); 
    tree->setInputCloud(cloud); 

    std::vector<pcl::PointIndices> clusterIndices; 
    pcl::EuclideanClusterExtraction<PointT> ec; 
    ec.setClusterTolerance(clusterTolerance); 
    ec.setMinClusterSize(minSize); 
    ec.setMaxClusterSize(maxSize); 
    ec.setSearchMethod(tree); 
    ec.setInputCloud(cloud); 
    ec.extract(clusterIndices); 

    for(pcl::PointIndices getIndices: clusterIndices)
    {
        typename pcl::PointCloud<PointT>::Ptr cloudCluster (new pcl::PointCloud<PointT>); 

        for(int index : getIndices.indices)
            cloudCluster->points.push_back (cloud->points[index]); 

        cloudCluster->width = cloudCluster->points.size(); 
        cloudCluster->height = 1; 
        cloudCluster->is_dense = true;

        clusters.push_back(cloudCluster); 

    }
    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "clustering took " << elapsedTime.count() << " milliseconds and found " << clusters.size() << " clusters" << std::endl;

    return clusters;
}

template<typename PointT>
void ProcessPointClouds<PointT>::clusterHelper(int indice, const std::vector<std::vector<float>> points, std::vector<int>& cluster,std::vector<bool>& processed, KdTree* tree, float distanceTol){
	processed[indice] = true; 
	cluster.push_back(indice); 

	std::vector<int> nearest = tree->search(points[indice], distanceTol); 

	for(int id : nearest){
		if(!processed[id])
			clusterHelper(id,points,cluster,processed,tree,distanceTol); 
	}
}
template<typename PointT>
std::vector<std::vector<int>> ProcessPointClouds<PointT>::euclideanCluster(const std::vector<std::vector<float>>& points, KdTree* tree, float distanceTol)
{

	// TODO: Fill out this function to return list of indices for each cluster
	std::vector<std::vector<int>> clusters;
    std::vector<bool> processed(points.size(),false); 

	int i =0; 
	while(i < points.size()){
		if(processed[i]){
			i++;
			continue; 
		}

		std::vector<int> cluster; 
		clusterHelper(i,points,cluster,processed,tree,distanceTol); 
		clusters.push_back(cluster); 
		i++; 


	}
	return clusters;

}

template<typename PointT>
std::vector<typename pcl::PointCloud<PointT>::Ptr> ProcessPointClouds<PointT>::Clustering3D(typename pcl::PointCloud<PointT>::Ptr cloud, float clusterTolerance, int minSize, int maxSize)
{
    // Time clustering process
    auto startTime = std::chrono::steady_clock::now();

    std::vector<typename pcl::PointCloud<PointT>::Ptr> clusters;
    std::vector<std::vector<float>> points; 
    // perform euclidean clustering to group detected obstacles
    KdTree* tree = new KdTree;
    for (int i =0; i<cloud->points.size(); i++){
        PointT point = cloud->points[i]; 
        tree->insert({point.x,point.y,point.z},i); 
    }
    
    std::vector<std::vector<int>> clusterInd = euclideanCluster(points,tree,clusterTolerance); 

    for(std::vector<int> cluster : clusterInd)
    {
        typename pcl::PointCloud<PointT>::Ptr cloudCluster (new pcl::PointCloud<PointT>); 

        for(int index : cluster)
            cloudCluster->points.push_back (cloud->points[index]); 

        cloudCluster->width = cloudCluster->points.size(); 
        cloudCluster->height = 1; 
        cloudCluster->is_dense = true;

        clusters.push_back(cloudCluster); 

    }
    auto endTime = std::chrono::steady_clock::now();
    auto elapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(endTime - startTime);
    std::cout << "clustering took " << elapsedTime.count() << " milliseconds and found " << clusters.size() << " clusters" << std::endl;

    return clusters;
}

template<typename PointT>
Box ProcessPointClouds<PointT>::BoundingBox(typename pcl::PointCloud<PointT>::Ptr cluster)
{

    // Find bounding box for one of the clusters
    PointT minPoint, maxPoint;
    pcl::getMinMax3D(*cluster, minPoint, maxPoint);

    Box box;
    box.x_min = minPoint.x;
    box.y_min = minPoint.y;
    box.z_min = minPoint.z;
    box.x_max = maxPoint.x;
    box.y_max = maxPoint.y;
    box.z_max = maxPoint.z;

    return box;
}


template<typename PointT>
void ProcessPointClouds<PointT>::savePcd(typename pcl::PointCloud<PointT>::Ptr cloud, std::string file)
{
    pcl::io::savePCDFileASCII (file, *cloud);
    std::cerr << "Saved " << cloud->points.size () << " data points to "+file << std::endl;
}


template<typename PointT>
typename pcl::PointCloud<PointT>::Ptr ProcessPointClouds<PointT>::loadPcd(std::string file)
{

    typename pcl::PointCloud<PointT>::Ptr cloud (new pcl::PointCloud<PointT>);

    if (pcl::io::loadPCDFile<PointT> (file, *cloud) == -1) //* load the file
    {
        PCL_ERROR ("Couldn't read file \n");
    }
    std::cerr << "Loaded " << cloud->points.size () << " data points from "+file << std::endl;

    return cloud;
}


template<typename PointT>
std::vector<boost::filesystem::path> ProcessPointClouds<PointT>::streamPcd(std::string dataPath)
{

    std::vector<boost::filesystem::path> paths(boost::filesystem::directory_iterator{dataPath}, boost::filesystem::directory_iterator{});

    // sort files in accending order so playback is chronological
    sort(paths.begin(), paths.end());

    return paths;

}