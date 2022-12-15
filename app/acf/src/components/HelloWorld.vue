<template>
  <h1>{{ msg }}</h1>

  <div class="card">
    <button type="button" @click="count++">count is {{ count }}</button>
    <p>
      Edit
      <code>components/HelloWorld.vue</code> to test HMRC
    </p>
  </div>

  <p class="read-the-docs">Click on the Vite and Vue logos to learn more</p>
</template>

<script>

export default {
  'name': "hw",
  props: {
    msg: String
  },
  data() {
    return {
      count: 0
    }
  }
}

</script>

<style scoped>
.read-the-docs {
  color: #888;
}
</style>
