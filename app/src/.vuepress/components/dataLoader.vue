<template>
<div>
  {{CurrentData}}
  Test
</div>
</template>

<script>
import pars_all from "../public/pars/pars_India.json";
import * as d3 from "d3";

const Locations = ["India", "Delhi", "Bihar", "Andhra_Pradesh"];

export default {
  name: "data_loader",
  data() {
    return {
      CurrentData: pars_all
    }
  },
  mounted() {
    console.log("mounted");
    this.fetchData(this.Location);
  },
  props: {
    Location: {
      type: String,
      default: "India"
    }
  },
  methods: {
    fetchData(loc) {
      console.log(loc)

    }
  }
}
</script>

<style scoped>

</style>