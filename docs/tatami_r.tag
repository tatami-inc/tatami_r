<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.12.0">
  <compound kind="file">
    <name>parallelize.hpp</name>
    <path>/github/workspace/include/tatami_r/</path>
    <filename>parallelize_8hpp.html</filename>
    <namespace>tatami_r</namespace>
  </compound>
  <compound kind="file">
    <name>sparse_matrix.hpp</name>
    <path>/github/workspace/include/tatami_r/</path>
    <filename>sparse__matrix_8hpp.html</filename>
    <namespace>tatami_r</namespace>
  </compound>
  <compound kind="file">
    <name>tatami_r.hpp</name>
    <path>/github/workspace/include/tatami_r/</path>
    <filename>tatami__r_8hpp.html</filename>
    <includes id="UnknownMatrix_8hpp" name="UnknownMatrix.hpp" local="yes" import="no" module="no" objc="no">UnknownMatrix.hpp</includes>
    <includes id="parallelize_8hpp" name="parallelize.hpp" local="yes" import="no" module="no" objc="no">parallelize.hpp</includes>
    <namespace>tatami_r</namespace>
  </compound>
  <compound kind="file">
    <name>UnknownMatrix.hpp</name>
    <path>/github/workspace/include/tatami_r/</path>
    <filename>UnknownMatrix_8hpp.html</filename>
    <includes id="parallelize_8hpp" name="parallelize.hpp" local="yes" import="no" module="no" objc="no">parallelize.hpp</includes>
    <class kind="struct">tatami_r::UnknownMatrixOptions</class>
    <class kind="class">tatami_r::UnknownMatrix</class>
    <namespace>tatami_r</namespace>
  </compound>
  <compound kind="class">
    <name>tatami_r::UnknownMatrix</name>
    <filename>classtatami__r_1_1UnknownMatrix.html</filename>
    <templarg>typename Value_</templarg>
    <templarg>typename Index_</templarg>
    <templarg>typename CachedValue_</templarg>
    <templarg>typename CachedIndex_</templarg>
    <base>tatami::Matrix&lt; Value_, Index_ &gt;</base>
    <member kind="function">
      <type></type>
      <name>UnknownMatrix</name>
      <anchorfile>classtatami__r_1_1UnknownMatrix.html</anchorfile>
      <anchor>a27b22842bfe04b6d5bf9a78fb4977102</anchor>
      <arglist>(Rcpp::RObject seed, const UnknownMatrixOptions &amp;opt)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UnknownMatrix</name>
      <anchorfile>classtatami__r_1_1UnknownMatrix.html</anchorfile>
      <anchor>abff9cf90f9cb485f30f60901f89e8ef8</anchor>
      <arglist>(Rcpp::RObject seed)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>tatami_r::UnknownMatrixOptions</name>
    <filename>structtatami__r_1_1UnknownMatrixOptions.html</filename>
    <member kind="variable">
      <type>std::optional&lt; std::size_t &gt;</type>
      <name>maximum_cache_size</name>
      <anchorfile>structtatami__r_1_1UnknownMatrixOptions.html</anchorfile>
      <anchor>ac231f00aa83f55652eb76299720617be</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>require_minimum_cache</name>
      <anchorfile>structtatami__r_1_1UnknownMatrixOptions.html</anchorfile>
      <anchor>a789336083111a71a161a58a1ad383c52</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>tatami_r</name>
    <filename>namespacetatami__r.html</filename>
    <class kind="class">tatami_r::UnknownMatrix</class>
    <class kind="struct">tatami_r::UnknownMatrixOptions</class>
    <member kind="function">
      <type>manticore::Executor &amp;</type>
      <name>executor</name>
      <anchorfile>namespacetatami__r.html</anchorfile>
      <anchor>a7cbde3ef2a02ed9a8d84ecd89697eb58</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_executor</name>
      <anchorfile>namespacetatami__r.html</anchorfile>
      <anchor>a8156d94510560edbd62605df13e033da</anchor>
      <arglist>(manticore::Executor *ptr)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parallelize</name>
      <anchorfile>namespacetatami__r.html</anchorfile>
      <anchor>a23b3149a67ca05913f08045bd1ef003c</anchor>
      <arglist>(const Function_ fun, const Index_ ntasks, int nthreads)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parse_SVT_SparseMatrix</name>
      <anchorfile>namespacetatami__r.html</anchorfile>
      <anchor>a1153096253856a330437d23e8f497d41</anchor>
      <arglist>(const Rcpp::RObject &amp;matrix, const Function_ fun)</arglist>
    </member>
  </compound>
  <compound kind="page">
    <name>md_parallel</name>
    <title>Enabling parallelization</title>
    <filename>md_parallel.html</filename>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>Read R objects via tatami</title>
    <filename>index.html</filename>
    <docanchor file="index.html">md__2github_2workspace_2README</docanchor>
  </compound>
</tagfile>
