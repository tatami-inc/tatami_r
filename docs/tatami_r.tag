<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile doxygen_version="1.10.0" doxygen_gitid="GIT-NOTFOUND">
  <compound kind="file">
    <name>parallelize.hpp</name>
    <path>/github/workspace/include/tatami_r/</path>
    <filename>parallelize_8hpp.html</filename>
    <namespace>tatami_r</namespace>
  </compound>
  <compound kind="file">
    <name>tatami_r.hpp</name>
    <path>/github/workspace/include/tatami_r/</path>
    <filename>tatami__r_8hpp.html</filename>
    <includes id="parallelize_8hpp" name="parallelize.hpp" local="yes" import="no" module="no" objc="no">parallelize.hpp</includes>
    <namespace>tatami_r</namespace>
  </compound>
  <compound kind="struct">
    <name>tatami_r::Options</name>
    <filename>structtatami__r_1_1Options.html</filename>
    <member kind="variable">
      <type>size_t</type>
      <name>maximum_cache_size</name>
      <anchorfile>structtatami__r_1_1Options.html</anchorfile>
      <anchor>a2692a265df2ac6b8160b5b5c6b68dcdc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>require_minimum_cache</name>
      <anchorfile>structtatami__r_1_1Options.html</anchorfile>
      <anchor>a5756e5c98133bfbeab0437024c13f1d4</anchor>
      <arglist></arglist>
    </member>
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
      <anchor>afb1de78326505497706b5278a66a1224</anchor>
      <arglist>(Rcpp::RObject seed, const Options &amp;opt)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>UnknownMatrix</name>
      <anchorfile>classtatami__r_1_1UnknownMatrix.html</anchorfile>
      <anchor>abff9cf90f9cb485f30f60901f89e8ef8</anchor>
      <arglist>(Rcpp::RObject seed)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>tatami_r</name>
    <filename>namespacetatami__r.html</filename>
    <class kind="struct">tatami_r::Options</class>
    <class kind="class">tatami_r::UnknownMatrix</class>
    <member kind="function">
      <type>manticore::Executor &amp;</type>
      <name>executor</name>
      <anchorfile>namespacetatami__r.html</anchorfile>
      <anchor>a7cbde3ef2a02ed9a8d84ecd89697eb58</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>parallelize</name>
      <anchorfile>namespacetatami__r.html</anchorfile>
      <anchor>aea54ff1c1c6a412b7c32f15814c0109e</anchor>
      <arglist>(Function_ fun, size_t njobs, size_t nthreads)</arglist>
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
