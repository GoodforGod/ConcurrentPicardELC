# Custom Picard ELC Implementation

A set of Java command line tools for manipulating high-throughput sequencing (HTS) data and formats.

All credits to Broad Institute [[ORIGINAL REPO](https://github.com/broadinstitute/picard)]

## What and Why

Main goal was to redesign **Estimated Library Complexity** from sequential to concurrent implementation, also SortingCollection was rewritten in _concurrent_ & _thread-safe_ way.   

## ELC Library Changes

There are 3 different concurrent implementation, the best and the most performance one **[STREAM]**, and the two least _code-clean_ and **may be** less _performance_.

Also **[EXECUTOR]** version contains all Predicates, abstract collection and etc, which are used in all Concurrent implementations.

P.S.
All _help_ and _support_ classes and predicates were grouped in **[EXECUTOR]** _intentionally_. 

_Estimated Library Complexity_
[[ORIGINAL](src/main/java/picard/sam/markduplicates/EstimateLibraryComplexity.java)]

_Custom Estimated Library Complexity_
[[STREAM](src/main/java/picard/sam/markduplicates/ConcurrentStreamedEstimateLibraryComplexity.java) | [POOL](src/main/java/picard/sam/markduplicates/ConcurrentPoolEstimateLibraryComplexity.java) | [EXECUTOR]( src/main/java/picard/sam/markduplicates/ConcurrentExecutorEstimateLibraryComplexity.java)]

### Sorting Collection

_Sorting Collection_
[[ORIGINAL](https://github.com/samtools/htsjdk/blob/master/src/main/java/htsjdk/samtools/util/SortingCollection.java)]

_Custom Sorting Collection_ 
[[CUSTOM](src/main/java/picard/sam/markduplicates/util/ConcurrentSortingCollection.java)]

### Other

Abstraction wrapper over PeekableIterator, used to read and filter sorted files in async mode.
[[QueueProducer](src/main/java/picard/sam/markduplicates/util/QueueProducer.java)]

## Performance Results

### First Test
| **Heap Size** | **Bam File Size** | **JDK**  | **Processor** |
| --- | --- | --- | --- | 
| *256mb* | *464mb* | *1.8.0_121 (64bit)* | *i3-4030u 1x2x2* |

#### Result

| **Implementation Name** | **Avegange Time (ms)** | **Amount of iteratins** | **MaxPairsInMemoty** |
| --- | --- | --- | --- |
| *Default*                 | *55541.009* | *5* | *169196* |
| *ConcurrentExecutorELC*   | *37775.565* | *5* | *169196* |
| *ConcurrentPoolELC*       | *39255.594* | *5* | *169196* |
| *ConcurrentStreamedELC*   | *40355.021* | *5* | *169196* |

* The best custom implementation is **47%** faster.

### Second Test
| **Heap Size** | **Bam File Size** | **JDK**  | **Processor** |
| --- | --- | --- | --- | 
| *4gb* | *16.9gb* | *1.8.0_121 (64bit)* | *i3-4030u 1x2x2* |

#### Result

| **Implementation Name** | **Avegange Time (ms)** | **Amount of iteratins** |
| --- | --- | --- |
| *Default*                 | *1972577.978* | *5* |
| *ConcurrentExecutorELC*   | *1253195.879* | *5* |

* Custom impelemntation is **57.4%** faster.

## Building Picard

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
