#Custom Picard ELC Implementation

A set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats.

All credits to Broad Institute [[ORIGINAL REPO](https://github.com/broadinstitute/picard)]

## What and Why

Main goal was to redesign **Estimated Library Complexity** from sequential to concurrent implementation, also SortingCollection was rewritten in _concurrent_ & _thread-safe_ way.   

### Differences

#### ELC

There are 3 different concurrent implementation, the best and the most performance one **[STREAM]**, and the two least _code-clean_ and **may be** less _performance_.

Also **[THE LEAST]** version contains all Predicates, abstract collection and etc, which are used in all Concurrent implementations.

P.S.
All _help_ and _support_ classes and predicates were grouped in **[EXECUTOR]** _intentionally_. 

**Estimated Library Complexity**
[[ORIGINAL]()]

**Custom Estimated Library Complexity** 
[[STREAM](src/main/java/picard/sam/markduplicates/ConcurrentStreamedEstimateLibraryComplexity.java) | [POOL](src/main/java/picard/sam/markduplicates/ConcurrentPoolEstimateLibraryComplexity.java) | [EXECUTOR]( src/main/java/picard/sam/markduplicates/ConcurrentExecutorEstimateLibraryComplexity.java)]

#### Sorting Collection

**Sorting Collection**
[[ORIGINAL](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/util/SortingCollection.java)]

**Custom Sorting Collection** 
[[CUSTOM](src/main/java/picard/sam/markduplicates/util/ConcurrentSortingCollection.java)]

#### Others

Abstraction wrapper over PeekableIterator, used to read and filter sorted files in async mode.
[[QueueProducer](src/main/java/picard/sam/markduplicates/util/QueueProducer.java)]

##Building Picard

* To build a fully-packaged, runnable Picard jar with all dependencies included, run:
```
    ./gradlew shadowJar
```

* The resulting jar will be in `build/libs`. To run it, the command is:
```
    java -jar build/libs/picard.jar
    
    or
    
    java -jar build/libs/picard-<VERSION>-all.jar 
```    

* To build a jar containing only Picard classes (without its dependencies), run:
```
    ./gradlew jar
```    
    
* To clean the build directory, run:
```
    ./gradlew clean
```

## Credits
[[ORIGINAL REPO](https://github.com/broadinstitute/picard)]